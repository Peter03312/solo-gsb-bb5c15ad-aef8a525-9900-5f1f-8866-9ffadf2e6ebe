"""双段（两份 chain 依次）迁移：引擎正确性与 HTTP 接口测试。"""

import os

import pytest
from fastapi.testclient import TestClient

from app.bed import BedRecord, parse_bed
from app.chain import Block, Chain, parse_chain
from app.main import app
from app.twoseg import TwoSegmentEngine

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _rec(start: int, end: int, chrom: str = "chrT", line_no: int = 1) -> BedRecord:
    return BedRecord(line_no, chrom, start, end, f"{chrom}\t{start}\t{end}", "")


def _chain1() -> Chain:
    # chrT + [0,100) -> chrM + [0,100)，单块
    return Chain(
        "1", "chrT", 200, "+", 0, 100,
        "chrM", 300, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )


def _chain2() -> Chain:
    # chrM + [0,100) -> chrF + [500,600)，单块
    return Chain(
        "2", "chrM", 300, "+", 0, 100,
        "chrF", 1000, "+", 500, 600, 1000,
        [Block(0, 100, 500, 600)], 1,
    )


def _engine(chains1=None, chains2=None, side1="target", side2="target"):
    return TwoSegmentEngine(
        chains1 if chains1 is not None else [_chain1()], side1,
        chains2 if chains2 is not None else [_chain2()], side2,
    )


# ---------- 基本映射与坐标回溯 ----------

def test_full_mapping_both_segments():
    row = _engine().convert(_rec(10, 40))
    assert row.status == "mapped"
    assert row.failed_stage == ""
    assert len(row.candidates) == 1
    path = row.candidates[0]
    assert (path.chain_id_1, path.chain_id_2) == ("1", "2")
    assert (path.covered1, path.covered2, path.total) == (30, 30, 30)
    frag = path.fragments[0]
    assert (frag.src_start, frag.src_end) == (10, 40)
    assert (frag.mid_start, frag.mid_end) == (10, 40)
    assert (frag.dest_start, frag.dest_end) == (510, 540)


def test_trace_back_through_minus_strand_segment1():
    # 第一段目标为负链：chrT + [0,100) -> chrM -（chrM 长 300，正向帧 [200,300)）
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 100,
        "chrM", 300, "-", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    # 第二段：chrM + [200,300) -> chrF + [0,100)
    c2 = Chain(
        "2", "chrM", 300, "+", 200, 300,
        "chrF", 1000, "+", 0, 100, 1000,
        [Block(200, 300, 0, 100)], 1,
    )
    eng = TwoSegmentEngine([c1], "target", [c2], "target")
    row = eng.convert(_rec(10, 40))
    assert row.status == "mapped"
    frag = row.candidates[0].fragments[0]
    # chrT [10,40) -> chrM 正向帧 [260,290) -> chrF [60,90)
    assert (frag.mid_start, frag.mid_end) == (260, 290)
    assert (frag.dest_start, frag.dest_end) == (60, 90)
    # 最终片段回溯到原 BED 源坐标（方向已镜像回来）
    assert (frag.src_start, frag.src_end) == (10, 40)
    assert frag.src_start < frag.src_end


def test_trace_back_partial_subrange_of_minus_strand_fragment():
    # 第二段只命中中间片段的子区间：回溯必须按映射方向取对应源子区间
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 100,
        "chrM", 300, "-", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    # 第二段只覆盖 chrM 正向帧 [270,300)（即源 chrT [0,30) 一段）
    c2 = Chain(
        "2", "chrM", 300, "+", 270, 300,
        "chrF", 1000, "+", 0, 30, 1000,
        [Block(270, 300, 0, 30)], 1,
    )
    eng = TwoSegmentEngine([c1], "target", [c2], "target")
    row = eng.convert(_rec(0, 100))
    assert row.status == "partial"
    path = row.candidates[0]
    assert (path.covered1, path.covered2, path.total) == (100, 30, 100)
    frag = path.fragments[0]
    # chrM [270,300) 对应源 chrT [0,30)（负链镜像回溯）
    assert (frag.src_start, frag.src_end) == (0, 30)
    assert (frag.mid_start, frag.mid_end) == (270, 300)
    assert (frag.dest_start, frag.dest_end) == (0, 30)
    assert "第二段" in row.reason and "70 bp" in row.reason  # 第二段丢失 70 bp


# ---------- 缺口与部分映射 ----------

def test_gap_in_segment1_keeps_uncovered_info():
    # 第一段 chrT [0,50) 覆盖、[50,100) 是缺口
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 50,
        "chrM", 300, "+", 0, 50, 1000,
        [Block(0, 50, 0, 50)], 1,
    )
    row = _engine(chains1=[c1]).convert(_rec(20, 80))
    assert row.status == "partial"
    path = row.candidates[0]
    assert (path.covered1, path.covered2, path.total) == (30, 30, 60)
    assert "30 bp" in row.reason
    # 残余片段仍按部分映射处理，绝不能标成 mapped
    assert row.status != "mapped"


def test_gap_in_segment2_marks_partial_not_mapped():
    # 第一段全覆盖；第二段只覆盖 chrM [0,60)
    c2 = Chain(
        "2", "chrM", 300, "+", 0, 60,
        "chrF", 1000, "+", 500, 560, 1000,
        [Block(0, 60, 500, 560)], 1,
    )
    row = _engine(chains2=[c2]).convert(_rec(10, 90))
    assert row.status == "partial"
    path = row.candidates[0]
    assert (path.covered1, path.covered2, path.total) == (80, 50, 80)
    assert (path.fragments[0].src_start, path.fragments[0].src_end) == (10, 60)
    assert "第二段" in row.reason and "30 bp" in row.reason


def test_gaps_in_both_segments_split_into_ordered_fragments():
    # 第一段：chrT [0,40)∪[60,100) -> chrM 同坐标；第二段：chrM [0,30)∪[70,100)
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 100,
        "chrM", 300, "+", 0, 100, 1000,
        [Block(0, 40, 0, 40), Block(60, 100, 60, 100)], 1,
    )
    c2 = Chain(
        "2", "chrM", 300, "+", 0, 100,
        "chrF", 1000, "+", 0, 60, 1000,
        [Block(0, 30, 0, 30), Block(70, 100, 30, 60)], 1,
    )
    row = _engine(chains1=[c1], chains2=[c2]).convert(_rec(0, 100))
    assert row.status == "partial"
    path = row.candidates[0]
    # 第一段丢 20（[40,60)），第二段再丢 20（[30,40)、[60,70)）
    assert (path.covered1, path.covered2, path.total) == (80, 60, 100)
    srcs = [(f.src_start, f.src_end) for f in path.fragments]
    assert srcs == [(0, 30), (70, 100)]          # 按源坐标有序
    dests = [(f.dest_start, f.dest_end) for f in path.fragments]
    assert dests == [(0, 30), (30, 60)]
    assert "第一段" in row.reason and "第二段" in row.reason


# ---------- 无映射与失败阶段 ----------

def test_side_parameters_apply_to_respective_segments():
    # side2="query"：中间坐标 chrM 位于第二份 chain 的 query 侧
    c2 = Chain(
        "2", "chrF", 1000, "+", 500, 600,
        "chrM", 300, "+", 0, 100, 1000,
        [Block(500, 600, 0, 100)], 1,
    )
    eng = TwoSegmentEngine([_chain1()], "target", [c2], "query")
    row = eng.convert(_rec(10, 40))
    assert row.status == "mapped"
    frag = row.candidates[0].fragments[0]
    assert (frag.dest_start, frag.dest_end) == (510, 540)
    assert row.candidates[0].dest_chrom == "chrF"

    # side1="query"：BED 位于第一份 chain 的 query 侧（chrM -> chrT -> chrF）
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 100,
        "chrM", 300, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    # 中间坐标是 chrT（chain1 的 target 侧），chain2 源侧必须是 chrT
    c2t = Chain(
        "2", "chrT", 200, "+", 0, 100,
        "chrF", 1000, "+", 500, 600, 1000,
        [Block(0, 100, 500, 600)], 1,
    )
    eng2 = TwoSegmentEngine([c1], "query", [c2t], "target")
    row2 = eng2.convert(_rec(10, 40, chrom="chrM"))
    assert row2.status == "mapped"
    # chrM [10,40) -> chrT [10,40) -> chrF [510,540)
    f2 = row2.candidates[0].fragments[0]
    assert (f2.src_start, f2.src_end) == (10, 40)
    assert (f2.mid_start, f2.mid_end) == (10, 40)
    assert (f2.dest_start, f2.dest_end) == (510, 540)


# ---------- 无映射与失败阶段 ----------

def test_unmapped_segment1_chrom_missing():
    row = _engine().convert(_rec(0, 10, chrom="chrX"))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment1"
    assert "第一段" in row.reason and "chrX" in row.reason
    assert row.candidates == []


def test_unmapped_segment1_gap():
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 50,
        "chrM", 300, "+", 0, 50, 1000,
        [Block(0, 50, 0, 50)], 1,
    )
    row = _engine(chains1=[c1]).convert(_rec(60, 90))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment1"
    assert "缺口" in row.reason or "对齐范围" in row.reason


def test_unmapped_segment2_chrom_missing():
    # 第二份 chain 的源侧没有 chrM —— 第二段无匹配源染色体，按该段无映射处理
    c2 = Chain(
        "2", "chrOther", 300, "+", 0, 100,
        "chrF", 1000, "+", 500, 600, 1000,
        [Block(0, 100, 500, 600)], 1,
    )
    row = _engine(chains2=[c2]).convert(_rec(10, 40))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment2"
    assert "第二段" in row.reason and "chrM" in row.reason
    assert "无一到达最终坐标" in row.reason


def test_unmapped_segment2_gap():
    c2 = Chain(
        "2", "chrM", 300, "+", 60, 100,
        "chrF", 1000, "+", 500, 540, 1000,
        [Block(60, 100, 500, 540)], 1,
    )
    row = _engine(chains2=[c2]).convert(_rec(0, 50))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment2"
    assert "第二段" in row.reason


# ---------- 多路径歧义：全部保留，不选择、不合并 ----------

def test_ambiguous_paths_all_kept():
    # 第一段两个 chain 都命中 chrT [0,100)
    c1b = Chain(
        "7", "chrT", 200, "+", 0, 100,
        "chrM2", 300, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    c2b = Chain(
        "8", "chrM2", 300, "+", 0, 100,
        "chrG", 1000, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    eng = TwoSegmentEngine(
        [_chain1(), c1b], "target", [_chain2(), c2b], "target"
    )
    row = eng.convert(_rec(10, 40))
    assert row.status == "ambiguous"
    assert len(row.candidates) == 2
    combos = {(p.chain_id_1, p.chain_id_2) for p in row.candidates}
    assert combos == {("1", "2"), ("7", "8")}
    assert [p.path_index for p in row.candidates] == [1, 2]
    assert "不自动选择或合并" in row.reason


def test_ambiguity_from_segment2_only():
    # 第一段唯一，第二段两个 chain 都命中中间区间
    c2b = Chain(
        "9", "chrM", 300, "+", 0, 100,
        "chrG", 1000, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    row = _engine(chains2=[_chain2(), c2b]).convert(_rec(10, 40))
    assert row.status == "ambiguous"
    assert len(row.candidates) == 2
    dests = sorted(p.dest_chrom for p in row.candidates)
    assert dests == ["chrF", "chrG"]


def test_dead_segment1_candidate_collapses_with_note():
    # 第一段两个候选，其中一个的中间染色体在第二段不存在 -> 只剩一条可达路径
    c1b = Chain(
        "7", "chrT", 200, "+", 0, 100,
        "chrNowhere", 300, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    row = _engine(chains1=[_chain1(), c1b]).convert(_rec(10, 40))
    assert row.status == "mapped"          # 唯一可达路径，不算歧义
    assert len(row.candidates) == 1
    assert "未形成可达路径" in row.reason


def test_same_chain2_shared_by_two_fragments_forms_one_path():
    # 第一段同一 chain 的两个片段都经同一第二段 chain：合并为一条路径
    c1 = Chain(
        "1", "chrT", 200, "+", 0, 100,
        "chrM", 300, "+", 0, 100, 1000,
        [Block(0, 40, 0, 40), Block(60, 100, 60, 100)], 1,
    )
    row = _engine(chains1=[c1]).convert(_rec(0, 100))
    assert row.status == "partial"
    assert len(row.candidates) == 1
    assert len(row.candidates[0].fragments) == 2


# ---------- 示例文件端到端 ----------

def test_example_files_end_to_end():
    bed = parse_bed(open(os.path.join(EXAMPLES, "example_twoseg.bed")).read())
    c1 = parse_chain(open(os.path.join(EXAMPLES, "demo.chain")).read())
    c2 = parse_chain(open(os.path.join(EXAMPLES, "demo2.chain")).read())
    assert bed.errors == [] and c1.errors == [] and c2.errors == []
    rows = TwoSegmentEngine(c1.chains, "target", c2.chains, "target").convert_all(
        bed.records
    )
    by_line = {r.line_no: r for r in rows}
    assert by_line[7].status == "mapped"
    f = by_line[7].candidates[0].fragments[0]
    assert (f.src_start, f.src_end, f.dest_start, f.dest_end) == (56, 176, 0, 120)

    assert by_line[8].status == "partial"
    assert by_line[8].candidates[0].covered1 == 126
    assert by_line[8].candidates[0].covered2 == 126

    assert by_line[9].status == "ambiguous"
    assert len(by_line[9].candidates) == 2

    assert by_line[10].status == "unmapped"
    assert by_line[10].failed_stage == "segment2"

    assert by_line[11].status == "unmapped"
    assert by_line[11].failed_stage == "segment1"


# ---------- HTTP 接口 ----------

def _post2(bed_name="example_twoseg.bed", side1="target", side2="target",
           chain1_name="demo.chain", chain2_name="demo2.chain",
           bed_bytes=None, chain1_bytes=None, chain2_bytes=None):
    def _payload(name, raw):
        if raw is not None:
            return (name, raw, "text/plain")
        return (name, open(os.path.join(EXAMPLES, name), "rb"), "text/plain")

    files = {
        "bed": _payload(bed_name, bed_bytes),
        "chain1": _payload(chain1_name, chain1_bytes),
        "chain2": _payload(chain2_name, chain2_bytes),
    }
    try:
        return client.post(
            "/api/liftover2", files=files,
            data={"side1": side1, "side2": side2},
        )
    finally:
        for f in files.values():
            f[1].close() if hasattr(f[1], "close") else None


def test_liftover2_preview_and_download():
    resp = _post2()
    assert resp.status_code == 200, resp.text
    data = resp.json()
    preview = data["preview"]
    assert preview["total"] == 5
    assert preview["counts"] == {
        "mapped": 1, "partial": 1, "ambiguous": 1, "unmapped": 2,
    }
    rows = {r["line_no"]: r for r in preview["rows"]}
    # 预览标明两段 chain ID、路径序号与每段覆盖情况
    amb = rows[9]
    assert [c["path_index"] for c in amb["candidates"]] == [1, 2]
    assert {c["chain_id_1"] for c in amb["candidates"]} == {"1", "2"}
    assert {c["chain_id_2"] for c in amb["candidates"]} == {"10", "11"}
    part = rows[8]["candidates"][0]
    assert (part["covered1"], part["covered2"], part["total"]) == (126, 126, 135)
    # 失败阶段
    assert rows[10]["failed_stage"] == "segment2"
    assert rows[11]["failed_stage"] == "segment1"

    token = data["download_token"]
    dl = client.get(f"/api/download/{token}")
    assert dl.status_code == 200
    body = dl.content.decode("utf-8-sig")
    header = body.splitlines()[0].split("\t")
    assert header == [
        "source_line", "status", "source_chrom", "source_start", "source_end",
        "path_index", "chain_id_seg1", "chain_id_seg2",
        "mid_chrom", "mid_strand", "dest_chrom", "dest_strand",
        "covered_seg1_bp", "covered_seg2_bp", "source_len",
        "fragment_index", "src_start", "src_end",
        "mid_start", "mid_end", "dest_start", "dest_end",
        "failed_stage", "reason",
    ]
    lines = body.splitlines()[1:]
    # 歧义行：两条路径各占一行，chain ID 成对出现
    amb_lines = [ln.split("\t") for ln in lines if "\tambiguous\t" in ln]
    assert len(amb_lines) == 2
    assert {(f[5], f[6], f[7]) for f in amb_lines} == {("1", "1", "10"), ("2", "2", "11")}
    # 无映射行保留源物理行号与失败阶段
    unm = {f[0]: f for f in (ln.split("\t") for ln in lines) if f[1] == "unmapped"}
    assert unm["10"][22] == "segment2"
    assert unm["11"][22] == "segment1"
    assert unm["10"][5] == "." and unm["11"][5] == "."
    # 部分映射行：每段覆盖情况
    part_line = next(ln.split("\t") for ln in lines if "\tpartial\t" in ln)
    assert (part_line[12], part_line[13], part_line[14]) == ("126", "126", "135")


def test_liftover2_bad_bed_rejected_no_download():
    resp = _post2(bed_bytes=b"chr1\t10\n")
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("第 1 行" in e for e in errors)
    assert "download_token" not in resp.json()


def test_liftover2_bad_chain1_rejected_with_line_number():
    bad_chain = b"chain 1000 chr1 100 + 0 10\n10\n"   # 头列数不足
    resp = _post2(chain1_bytes=bad_chain)
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("chain1" in e and "第 1 行" in e for e in errors)


def test_liftover2_bad_chain2_rejected_with_line_number():
    bad_chain = b"chain 1000 chrA 500 + 0 20 chrM 600 + 0 20 1\n15\n"  # 终止块未到末端
    resp = _post2(chain2_bytes=bad_chain)
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("chain2" in e and "第 2 行" in e for e in errors)


def test_liftover2_invalid_side_rejected():
    resp = _post2(side2="upside")
    assert resp.status_code == 400
    assert "side2" in resp.json()["errors"][0]


def test_liftover2_job_not_adjudicable():
    token = _post2().json()["download_token"]
    assert client.get(f"/api/result/{token}/decisions").status_code == 400
    resp = client.post(
        f"/api/result/{token}/decisions",
        json={"version": 0, "source_line": 9, "candidate_index": 1},
    )
    assert resp.status_code == 400
    assert client.get(f"/api/download/{token}/adjudicated").status_code == 400
    # 原始下载不受影响
    assert client.get(f"/api/download/{token}").status_code == 200


def test_single_segment_still_adjudicable():
    # 单段流程不受影响（兼容性回归）
    with open(os.path.join(EXAMPLES, "example_target.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover",
            files={"bed": ("b.bed", bf, "text/plain"),
                   "chain": ("c.chain", cf, "text/plain")},
            data={"side": "target"},
        )
    assert resp.status_code == 200
    token = resp.json()["download_token"]
    assert client.get(f"/api/result/{token}/decisions").status_code == 200
    assert client.get(f"/api/download/{token}/adjudicated").status_code == 200
