# 离线基因组区间迁移 Web 应用（BED Liftover）

供实验室**内网/离线环境**使用的基因组坐标区间迁移工具。用户在浏览器中
上传 BED 文件与对应的 UCSC `chain` 文件，指定 BED 区间属于链的
**target 侧**还是**query 侧**，先在页面上逐条预览转换状态，确认后再
下载 TSV 结果。

支持两种转换模式：

- **单段**：一份 chain，BED 区间从链的一侧迁到另一侧；
- **双段**：同一 BED 依次经过**两份 chain** 迁移（如 hg19→hg38→T2T），
  第二段只处理第一段实际映射到的片段，最终片段按映射方向精确回溯到
  原 BED 源坐标，并标明两段 chain ID、路径序号与每段覆盖情况。

- 后端：Python 3.11+ / FastAPI / Uvicorn
- 前端：原生 HTML + CSS + JavaScript（无任何 CDN / 外网依赖）
- 运行时数据不落盘：上传内容与转换结果只保存在进程内存中；校验失败时
  **不会**生成任何可下载文件

---

## 目录结构

```
<交付外层目录>/
├── Dockerfile            # 方式 A 构建入口：以外层目录为上下文（COPY liftover-web/...）
├── .dockerignore         # 外层构建上下文排除规则
└── liftover-web/
    ├── app/                  # 后端
    │   ├── bed.py            # BED 解析（零基左闭右开，错误带行号）
    │   ├── chain.py          # UCSC chain 解析（严格校验，错误带行号）
    │   ├── engine.py         # 单段区间转换引擎（缺口拆分 / 歧义 / 负链）
    │   ├── twoseg.py         # 双段迁移引擎（路径组合 / 坐标回溯 / 失败阶段）
    │   ├── report.py         # 单段原始 / 裁决版与双段 TSV 导出、预览 JSON
    │   └── main.py           # FastAPI 入口（单段 + 双段 + 歧义人工裁决接口）
    ├── static/               # 前端（纯本地资源）
    ├── examples/             # 可直接试用的 BED / chain 示例（含双段示例）
    ├── tests/                # pytest 测试（68 个用例）
    ├── scripts/init.sh       # 本地开发依赖初始化
    ├── Dockerfile            # 方式 B 构建入口：以本目录为上下文
    ├── docker-compose.yml    # 两种启动方式共用
    ├── .env.example
    ├── requirements.txt      # 运行依赖
    └── requirements-dev.txt  # 运行 + 测试依赖
```

---

## 一、Docker / Compose 启动（推荐，完全离线）

镜像在构建时安装全部 Python 依赖，运行阶段不需要访问外网。下面两种
入口构建的是**同一个镜像** `liftover-web:1.0.0`：应用代码与静态资源
全部打进镜像，缺少任一构建必需文件（`requirements.txt`、`app/`、
`static/`）时构建直接失败，不会启动缺少资源却显示健康的容器。

### 方式 A：在外层目录一键启动

在交付外层目录直接指定 `liftover-web/docker-compose.yml`，并以外层
目录为 Compose 项目目录：

```bash
# 构建并后台启动
docker compose --project-directory . -f liftover-web/docker-compose.yml up -d --build

# 查看状态 / 日志
docker compose --project-directory . -f liftover-web/docker-compose.yml ps
docker compose --project-directory . -f liftover-web/docker-compose.yml logs -f

# 停止
docker compose --project-directory . -f liftover-web/docker-compose.yml down
```

此时构建上下文是外层目录，使用外层目录下的 `Dockerfile`（COPY 源路径
带 `liftover-web/` 前缀）。

### 方式 B：进入 liftover-web 目录启动

```bash
cd liftover-web

# 构建并后台启动
docker compose up -d --build

# 查看状态 / 日志
docker compose ps
docker compose logs -f

# 停止
docker compose down
```

此时构建上下文是 `liftover-web/` 目录，使用该目录下的 `Dockerfile`。
两份 Dockerfile 仅 COPY 源路径前缀不同，构建产物一致。

访问地址（两种方式相同）：

- **http://localhost:8080** （应用页面）
- 健康检查：`http://localhost:8080/healthz`

### 端口与配置

复制环境变量模板并按需修改（方式 B 在 `liftover-web/` 内执行；方式 A
可将 `.env` 放在外层目录，或直接在命令行注入变量）：

```bash
cp .env.example .env
```

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `LIFTOVER_PORT` | `8080` | 宿主机映射端口（容器内固定监听 8080） |
| `LIFTOVER_MAX_UPLOAD_MB` | `50` | 单个上传文件大小上限（MB） |
| `LIFTOVER_JOB_TTL_SECONDS` | `3600` | 转换结果在内存中的保留时长（秒） |
| `LIFTOVER_MAX_JOBS` | `100` | 内存中同时保留的结果任务数量上限 |

改端口示例（两种入口写法都列出来）：

```bash
# 方式 A（外层目录）
LIFTOVER_PORT=9000 docker compose --project-directory . -f liftover-web/docker-compose.yml up -d

# 方式 B（liftover-web/ 内）
LIFTOVER_PORT=9000 docker compose up -d
```

随后访问 `http://localhost:9000`。

不用 Compose 时：

```bash
# 在 liftover-web/ 内构建
docker build -t liftover-web:1.0.0 .
# 或在外层目录构建（等价）
docker build -f Dockerfile -t liftover-web:1.0.0 .

docker run -d --name liftover-web -p 8080:8080 \
  -e LIFTOVER_MAX_UPLOAD_MB=100 \
  liftover-web:1.0.0
```

---

## 二、本地 Python 直接启动（开发 / 调试）

需要 Python 3.11+。

```bash
cd liftover-web

# 1) 初始化：创建虚拟环境并安装依赖（含测试依赖）
bash scripts/init.sh

# 2) 启动
.venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8080
```

访问 **http://localhost:8080**。

运行测试：

```bash
.venv/bin/python -m pytest -q
```

---

## 三、使用方法

### 单段转换

1. 转换模式选择「单段（一份 chain）」（默认）。
2. 选择 BED 文件（至少前 3 列 `chrom start end`，坐标为**零基左闭右开**，
   支持第 4 列起的额外列但不参与转换；允许 `#`、`browser`、`track` 行）。
3. 选择对应的 UCSC `chain` 文件。
4. 选择 BED 区间位于 chain 的哪一侧：
   - **target 侧**：target → query；
   - **query 侧**：query → target。
5. 点击「校验并转换」。页面显示每条源区间的状态。
6. 若存在**歧义**行，可先在页面上按 BED 物理行号逐行**人工裁决**
   （见下一节），也可以直接下载。
7. 点击「下载原始结果 TSV」获取未改动的原始结果；对歧义行做过裁决后，
   可点击「下载裁决版 TSV」获取只保留所选候选的裁决版。

### 双段转换（两份 chain 依次迁移）

1. 转换模式选择「双段（两份 chain 依次迁移）」。
2. 选择 BED 文件与 **chain 文件 ①（第一段）**，并指定 BED 区间属于
   chain1 的哪一侧（`side1`）。
3. 选择 **chain 文件 ②（第二段）**，并指定第一段映射到的**中间坐标**
   属于 chain2 的哪一侧（`side2`）。例如 hg19→hg38（chain1，BED 在其
   target 侧）再接 hg38→T2T（chain2，中间坐标 hg38 在其 target 侧），
   两段都选 target 侧。
4. 点击「校验并转换」。第二段**只处理第一段实际映射到的片段**；每个
   最终片段都按映射方向（含负链镜像）回溯到原 BED 源坐标。
5. 预览中每条路径标明**两段 chain ID、路径序号与每段覆盖（bp）**，
   片段按「源 → 中间 → 最终」坐标展开；点击「下载原始结果 TSV」获取
   双段结果。

双段结果的约定：

- **多条可达路径全部保留为歧义**（`ambiguous`），不自动选择、不合并；
  路径 = 第一段某 chain 命中 × 第二段某 chain 命中的组合，同一组合跨
  多个中间片段时合并为一条路径；
- **任一段的缺口都会保留源区间未覆盖信息**：`covered_seg1_bp` /
  `covered_seg2_bp` 分别记录第一段覆盖的源 bp 与最终到达的源 bp，
  残余片段绝不标记为完全映射（`mapped`）；
- **无最终映射的源行**保留原物理行号与失败阶段（`failed_stage` 为
  `segment1` / `segment2`）；第二段无匹配源染色体同样按该段无映射处理；
- 双段结果不提供裁决版下载（裁决接口仅适用于单段结果，对双段 token
  调用会返回 400）。

### 歧义区间人工裁决（单段）

歧义区间（`ambiguous`）不能自动进入下游分析，必须经人工裁决。
裁决**只作用于这一次转换结果**，不会改动原始转换：

- 在预览表的每个歧义行，候选按候选序号（1、2、…，与 TSV 的
  `candidate_index` 一致）列出，点候选前的单选框即提交裁决；
  点「撤销选择」可恢复为未裁决（重新输出全部候选）；
- 页面上方实时显示「已选/未裁决行数」与**结果版本号**；
- 同一个结果被连续修改时版本号递增；若两人/两个标签页同时拿着同一
  版本提交，只有一次成功，另一次收到 409 并自动刷新当前选择后重试；
- 原始下载（`GET /api/download/{token}`）的内容在裁决前后**字节不变**，
  行的原始 `status`（如 `ambiguous`/`partial`）也始终保留——裁决不会把
  部分映射改称完全映射。

### 状态定义

| 状态 | 含义 |
| --- | --- |
| 完全映射 `mapped` | 区间全部落在对齐块内（双段：全部源 bp 经两段到达最终坐标） |
| 部分映射 `partial` | 一部分落在链内缺口/未对齐区域；命中部分被拆成**按源坐标有序的片段**，未映射部分丢弃（双段：任一段的缺口都会按源区间保留未覆盖 bp 数，残余片段不视为完全映射） |
| 歧义 `ambiguous` | 同一源区间命中**多个 chain / 多个目标位置**（双段：存在多条可达路径）；系统列出全部候选，**绝不自动挑选或合并**，单段结果需人工逐行裁决后再用于下游分析 |
| 无映射 `unmapped` | 没有任何碱基可映射；结果中保留源行号、源区间与原因（双段另以 `failed_stage` 标明失败于第一段还是第二段） |

### 坐标与负链

- 输入与输出一律是**零基、左闭右开**（`start < end`）。
- chain 负链侧的存储坐标按 UCSC 语义（以染色体长度为枢轴）归一化后
  再做区间求交，因此**负链转换后的坐标仍然是 start < end 的零基左闭
  右开区间**；结果的 `dest_strand` 列标明该候选在目标坐标系中的方向。

### 下载文件格式（TSV）

UTF-8（带 BOM，Excel 友好）。

**原始结果**（`下载原始结果 TSV` / `GET /api/download/{token}`）列：

```
source_line, status, source_chrom, source_start, source_end,
candidate_index, chain_id, dest_chrom, dest_strand,
fragment_index, src_start, src_end, dest_start, dest_end, reason
```

- `source_line` 是源文件**物理行号**（从 1 开始）；
- 跨缺口区间每个片段占一行；歧义区间的每个候选、每个片段各占一行；
- `unmapped` 行保留源行号与原因，候选/片段列写 `.`。

**裁决版**（`下载裁决版 TSV` / `GET /api/download/{token}/adjudicated`）
在原始 15 列末尾追加一列 `decision`，取值：

| decision | 出现在 | 输出内容 |
| --- | --- | --- |
| `已裁决` | 已选定候选的歧义行 | **仅输出被选候选的全部片段**（`candidate_index` 沿用原序号） |
| `未裁决` | 尚未选择的歧义行 | 仍输出**全部候选**的全部片段 |
| `.` | mapped / partial / unmapped 行 | 与原始结果对应行**逐列一致** |

裁决版的 `status` 列与原始结果完全一致：歧义行即使选了候选仍写作
`ambiguous`，`partial` 也不会因为人工选择而变成 `mapped`。
文件名形如 `liftover_target_to_query_<时间戳>.adjudicated.v<版本>.tsv`。

**双段结果**（`POST /api/liftover2` 后 `GET /api/download/{token}`）列：

```
source_line, status, source_chrom, source_start, source_end,
path_index, chain_id_seg1, chain_id_seg2,
mid_chrom, mid_strand, dest_chrom, dest_strand,
covered_seg1_bp, covered_seg2_bp, source_len,
fragment_index, src_start, src_end, mid_start, mid_end, dest_start, dest_end,
failed_stage, reason
```

- `path_index` 是该源行内 1 基路径序号；`chain_id_seg1` / `chain_id_seg2`
  是两段各自命中的 chain ID；
- `covered_seg1_bp` 为第一段覆盖的源 bp，`covered_seg2_bp` 为经两段后
  到达最终坐标的源 bp，`source_len` 为源区间长度；
- `src_start/src_end` 是回溯到原 BED 源坐标的片段范围，`mid_*` 为中间
  坐标，`dest_*` 为最终坐标（均为零基左闭右开）；
- `unmapped` 行路径/片段列写 `.`，`failed_stage` 为 `segment1` 或
  `segment2`，其余行该列写 `.`。

### 失败即阻止导出

- BED 坏行（列数不足、坐标非整数、`start >= end` 等）或 chain 格式错误
  （头列数/字段非法、块越界、缺少终止行、gap 越界等）都会在页面上
  **带行号列出全部问题**，接口返回 400，且**不会生成或保留任何可下载
  文件**；修正后需重新提交。双段模式下 BED 与两份 chain **任一文件**
  出错都会拒绝整次转换，错误信息以 `chain1 文件：` / `chain2 文件：`
  前缀标明来源。

---

## 四、HTTP API 调用示例

页面的全部功能都通过以下 HTTP 接口完成（端口以 `LIFTOVER_PORT` 为准，
默认 8080；示例中的文件路径按方式 A 的外层目录给出，方式 B 下去掉
`liftover-web/` 前缀即可）。

### 健康检查

```bash
curl http://localhost:8080/healthz
# {"status":"ok"}
```

### 上传并转换（预览）

`POST /api/liftover`，`multipart/form-data`，字段：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `bed` | 文件 | BED 文件（零基左闭右开） |
| `chain` | 文件 | UCSC chain 文件 |
| `side` | 文本 | `target`（target → query）或 `query`（query → target） |

```bash
curl -s -X POST http://localhost:8080/api/liftover \
  -F "bed=@liftover-web/examples/example_target.bed" \
  -F "chain=@liftover-web/examples/demo.chain" \
  -F "side=target"
```

成功（200）返回预览与下载令牌：

```json
{
  "download_token": "……",
  "download_name": "liftover_target_to_query_20261005-120000.tsv",
  "preview": {
    "total": 6,
    "counts": {"mapped": 1, "partial": 1, "ambiguous": 1, "unmapped": 3},
    "rows": []
  }
}
```

输入校验失败（400）只返回带行号的错误列表，**不会**产生任何可下载
文件：

```json
{"errors": ["BED 第 3 行：……", "……"]}
```

### 双段转换（两份 chain 依次迁移）

`POST /api/liftover2`，`multipart/form-data`，字段：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `bed` | 文件 | BED 文件（零基左闭右开） |
| `chain1` | 文件 | 第一份 UCSC chain 文件（第一段） |
| `side1` | 文本 | BED 区间在 chain1 中的源侧：`target` 或 `query` |
| `chain2` | 文件 | 第二份 UCSC chain 文件（第二段） |
| `side2` | 文本 | 第一段映射到的中间坐标在 chain2 中的源侧：`target` 或 `query` |

```bash
curl -s -X POST http://localhost:8080/api/liftover2 \
  -F "bed=@liftover-web/examples/example_twoseg.bed" \
  -F "chain1=@liftover-web/examples/demo.chain"  -F "side1=target" \
  -F "chain2=@liftover-web/examples/demo2.chain" -F "side2=target"
```

成功（200）返回结构与单段一致（`download_token` / `download_name` /
`preview`），预览的每个候选是一条**路径**，标明两段 chain ID、路径
序号与每段覆盖（`covered1` / `covered2` / `total`），片段含
`src_*`（回溯后的源坐标）、`mid_*`（中间坐标）、`dest_*`（最终坐标）；
无最终映射的行带 `failed_stage`（`segment1` / `segment2`）：

```json
{
  "download_token": "……",
  "download_name": "liftover_twoseg_20261005-120000.tsv",
  "preview": {
    "total": 5,
    "counts": {"mapped": 1, "partial": 1, "ambiguous": 1, "unmapped": 2},
    "rows": [
      {"line_no": 9, "status": "ambiguous", "failed_stage": "",
       "candidates": [
         {"path_index": 1, "chain_id_1": "1", "chain_id_2": "10",
          "mid_chrom": "chrA", "mid_strand": "+",
          "dest_chrom": "chrM", "dest_strand": "+",
          "covered1": 40, "covered2": 40, "total": 40,
          "fragments": [{"src_start": 320, "src_end": 360,
                          "mid_start": 364, "mid_end": 404,
                          "dest_start": 264, "dest_end": 304}]},
         {"path_index": 2, "chain_id_1": "2", "chain_id_2": "11", "……": "……"}
       ]}
    ]
  }
}
```

- 任一文件（BED / chain1 / chain2）格式错误都返回 400 并带行号
  （chain 错误以 `chain1 文件：` / `chain2 文件：` 前缀标明来源），
  整次转换被拒绝，不产生下载令牌；
- 返回的 `download_token` 可照常用于 `GET /api/download/{token}`
  下载双段 TSV（24 列，见「下载文件格式」一节）；
- 双段结果的多条可达路径全部保留、不做人工裁决：对该 token 调用
  裁决接口（`GET/POST /api/result/{token}/decisions`、
  `GET /api/download/{token}/adjudicated`）返回 400。

### 下载结果 TSV

`GET /api/download/{download_token}`。token 来自上一步响应，只在转换
完全成功后才存在，且结果只保存在服务进程内存中：

```bash
TOKEN=<上一步返回的 download_token>
curl -OJ "http://localhost:8080/api/download/${TOKEN}"
```

- 成功（200）：`text/tab-separated-values`（UTF-8 带 BOM），
  `Content-Disposition` 携带上一步返回的 `download_name` 文件名；
- token 不存在或已过期（404）：
  `{"errors": ["结果不存在或已过期，请重新提交转换"]}`。

### 查看当前裁决选择

`GET /api/result/{download_token}/decisions`。返回某次结果的版本号与
全部歧义行的当前选择（`candidate_index` 为 `null` 表示未裁决）：

```bash
curl -s "http://localhost:8080/api/result/${TOKEN}/decisions"
```

```json
{
  "download_token": "……",
  "version": 0,
  "decisions": [
    {"source_line": 11, "candidate_index": null, "candidate_count": 2}
  ]
}
```

### 提交人工裁决（选定候选 / 撤销）

`POST /api/result/{download_token}/decisions`，`application/json`：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `version` | 整数 | 客户端所基于的**当前结果版本**（首次为 0） |
| `source_line` | 整数 | BED 源文件**物理行号**，该行必须是歧义行 |
| `candidate_index` | 整数或 `null` | 选定的 1 基候选序号；`null` 表示撤销选择 |

```bash
# 为第 11 行选定候选 2（基于版本 0）
curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/decisions" \
  -H 'Content-Type: application/json' \
  -d '{"version":0,"source_line":11,"candidate_index":2}'
```

成功修改（200）后版本递增：

```json
{"download_token":"……","version":1,"modified":true,
 "source_line":11,"candidate_index":2}
```

提交与当前选择相同（无实际变化）时同样返回 200，但 `modified` 为
`false`、`version` 不变。随后可在版本 1 上改选候选 1 或撤销：

```bash
# 撤销第 11 行的选择（基于版本 1）
curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/decisions" \
  -H 'Content-Type: application/json' \
  -d '{"version":1,"source_line":11,"candidate_index":null}'
```

明确拒绝的情形：

| 场景 | 状态码 | 说明 |
| --- | --- | --- |
| token 不存在 / 结果已过期 | 404 | 需重新提交转换 |
| `version` 与服务端当前版本不一致 | **409** | 相同版本的并发修改只有一次成功；响应体带 `current_version`，查看最新选择后用新版本重试 |
| `source_line` 不存在于本次结果 | 400 | 错误信息带行号 |
| 该行不是歧义行（mapped/partial/unmapped） | 400 | 只有 `ambiguous` 行允许裁决 |
| `candidate_index` 越界（0、负数或超过候选数） | 400 | 错误信息标明允许范围 1..N |
| 请求体不是 JSON 对象 / 缺字段 / 字段类型错误 | 400 | 如 `version` 传字符串、`candidate_index` 传 `true` |

### 下载裁决版 TSV

`GET /api/download/{download_token}/adjudicated`，无需额外参数，
始终按服务端当前裁决状态即时生成：

```bash
curl -OJ "http://localhost:8080/api/download/${TOKEN}/adjudicated"
# 文件名：liftover_target_to_query_<时间戳>.adjudicated.v<版本>.tsv
```

- 已选歧义行只输出被选候选的全部片段（`decision=已裁决`）；
- 未选歧义行输出全部候选并标明 `decision=未裁决`；
- 其他状态行照原结果输出（`decision=.`），`status` 列保持原值；
- 与原始下载互不影响：`GET /api/download/{token}` 仍是裁决前的内容。

---

## 五、快速试用

单段：

```bash
# 页面上分别选择：
#   模式:  单段（一份 chain）
#   BED:   examples/example_target.bed
#   chain: examples/demo.chain
#   侧别:  target 侧
```

该示例覆盖全部四种状态（完全映射 / 跨缺口拆分的部分映射 / 歧义 /
无映射）。`examples/example_query.bed` 演示 query → target，
包含负链 query 的反向转换；`examples/example_bad.bed` 演示错误定位
与导出拦截。

歧义行（示例 BED 第 11 行 `chr1 320 360`，候选 1=chrA、候选 2=chrB）
可直接在预览页裁决，或用上面的 `POST /api/result/{token}/decisions`
接口按物理行号选定 / 撤销，再下载 `…/adjudicated` 裁决版 TSV。

双段：

```bash
# 页面上分别选择：
#   模式:   双段（两份 chain 依次迁移）
#   BED:    examples/example_twoseg.bed
#   chain①: examples/demo.chain  （target 侧）
#   chain②: examples/demo2.chain （target 侧）
```

该示例覆盖双段的全部情形：两段完全映射、第一段缺口导致的部分映射、
两条可达路径并存的歧义（`chain 1→10` 与 `chain 2→11`，含负链第二段）、
失败于第二段（中间片段不在 chain2 覆盖内）与失败于第一段（源染色体
不在 chain1 源侧）的无映射行。

---

## 六、离线部署说明

- 前端资源全部由本服务在 `/static/` 下提供，页面不含任何外部字体、
  脚本或样式表；
- Docker 构建需在有网络的机器上 `docker build`（或导出镜像），之后把
  镜像 `docker save | gzip > liftover-web.tgz` 拷贝到内网，
  `docker load` 后即可运行；
- 容器以非 root 用户运行，只暴露 8080 端口，不挂载任何数据卷；
- 上传文件与结果仅存内存，进程重启即清除，适合处理受限数据。
