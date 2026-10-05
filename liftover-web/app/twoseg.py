"""双段（两份 chain 依次）区间迁移引擎。

第一段把 BED 源区间映射到中间坐标；第二段**只处理第一段实际映射到的
片段**，把中间片段映射到最终坐标。每个最终片段都按映射方向（含负链
镜像）精确回溯到原 BED 源坐标，输出仍为零基左闭右开。

约定：

- 任一段的缺口都会丢掉对应碱基：结果按源区间保留未覆盖信息
  （covered_seg1 / covered_seg2 与原因文本），残余片段绝不标记为
  完全映射；
- 多条可达路径全部保留为歧义（ambiguous），不自动选择、不合并；
- 无最终映射的源行保留原物理行号与失败阶段（``failed_stage`` 为
  ``segment1`` / ``segment2``）；第二段无匹配源染色体同样按该段
  无映射处理。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .bed import BedRecord
from .chain import Chain
from .engine import _SegmentIndex


@dataclass(frozen=True)
class FinalFragment:
    """双段迁移的最终片段，坐标均为零基左闭右开。"""

    src_start: int               # 原 BED 源坐标（回溯后）
    src_end: int
    mid_start: int               # 第一段目标（中间）坐标
    mid_end: int
    dest_start: int              # 第二段目标（最终）坐标
    dest_end: int


@dataclass
class PathCandidate:
    """一条可达路径：第一段某个 chain 命中 × 第二段某个 chain 命中。"""

    path_index: int              # 1 基，在该源行内的路径序号
    chain_id_1: str
    chain_id_2: str
    mid_chrom: str
    mid_strand: str              # 第一段目标侧方向
    dest_chrom: str
    dest_strand: str             # 第二段目标侧方向
    fragments: list[FinalFragment]
    covered1: int                # 第一段覆盖的源 bp
    covered2: int                # 经两段后到达最终坐标的源 bp
    total: int                   # 源区间长度


@dataclass
class TwoSegRowResult:
    line_no: int
    src_chrom: str
    src_start: int
    src_end: int
    raw: str
    status: str                  # mapped / partial / ambiguous / unmapped
    reason: str = ""
    failed_stage: str = ""       # unmapped 时：segment1 / segment2
    candidates: list[PathCandidate] = field(default_factory=list)


class TwoSegmentEngine:
    """按顺序串联两份 chain 的迁移引擎。"""

    def __init__(
        self,
        chains1: list[Chain],
        side1: str,
        chains2: list[Chain],
        side2: str,
    ) -> None:
        self.seg1 = _SegmentIndex(chains1, side1)
        self.seg2 = _SegmentIndex(chains2, side2)
        self.side1 = self.seg1.side
        self.side2 = self.seg2.side

    def convert(self, rec: BedRecord) -> TwoSegRowResult:
        row = TwoSegRowResult(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            raw=rec.raw,
            status="unmapped",
        )
        total = rec.end - rec.start

        hits1, chrom1_known = self.seg1.candidates_for(
            rec.chrom, rec.start, rec.end
        )
        if not hits1:
            row.failed_stage = "segment1"
            if not chrom1_known:
                row.reason = (
                    f"第一段：染色体 {rec.chrom} 在第一份 chain 文件的源侧不存在"
                )
            else:
                row.reason = (
                    f"第一段：[{rec.start},{rec.end}) 完全落在第一份 chain 的"
                    f"对齐范围之外（链内缺口或未覆盖区域），没有碱基进入第二段"
                )
            return row

        # 第二段只处理第一段实际映射到的片段；路径按 (段1 chain, 段2 chain)
        # 的 chain 对象组合分组，同一组合跨多个中间片段的片段合并到同一路径。
        paths: list[PathCandidate] = []
        path_by_key: dict[tuple[int, int], PathCandidate] = {}
        seg2_missing_chroms: set[str] = set()
        seg2_gap_frags = 0
        reached_hits1: set[int] = set()

        for h1 in hits1:
            for tf in h1.fragments:
                hits2, chrom2_known = self.seg2.candidates_for(
                    h1.pc.dest_chrom, tf.dest_start, tf.dest_end
                )
                if not hits2:
                    if not chrom2_known:
                        seg2_missing_chroms.add(h1.pc.dest_chrom)
                    else:
                        seg2_gap_frags += 1
                    continue
                reached_hits1.add(id(h1))
                for h2 in hits2:
                    key = (id(h1.pc), id(h2.pc))
                    path = path_by_key.get(key)
                    if path is None:
                        path = PathCandidate(
                            path_index=0,        # 收集完毕后统一编号
                            chain_id_1=h1.pc.chain_id,
                            chain_id_2=h2.pc.chain_id,
                            mid_chrom=h1.pc.dest_chrom,
                            mid_strand=h1.pc.dest_strand,
                            dest_chrom=h2.pc.dest_chrom,
                            dest_strand=h2.pc.dest_strand,
                            fragments=[],
                            covered1=h1.covered,
                            covered2=0,
                            total=total,
                        )
                        path_by_key[key] = path
                        paths.append(path)
                    for f2 in h2.fragments:
                        # 把第二段片段的中间源区间回溯到原 BED 源坐标
                        s_lo, s_hi = tf.invert_dest(f2.src_start, f2.src_end)
                        path.fragments.append(
                            FinalFragment(
                                src_start=s_lo,
                                src_end=s_hi,
                                mid_start=f2.src_start,
                                mid_end=f2.src_end,
                                dest_start=f2.dest_start,
                                dest_end=f2.dest_end,
                            )
                        )

        if not paths:
            row.failed_stage = "segment2"
            bits: list[str] = []
            if seg2_missing_chroms:
                bits.append(
                    "中间染色体 "
                    + "、".join(sorted(seg2_missing_chroms))
                    + " 在第二份 chain 文件的源侧不存在"
                )
            if seg2_gap_frags:
                bits.append(
                    f"{seg2_gap_frags} 个中间片段落在第二份 chain 的对齐范围之外"
                    f"（链内缺口或未覆盖区域）"
                )
            row.reason = (
                "第二段：" + "；".join(bits)
                + "，第一段已映射的碱基无一到达最终坐标"
            )
            return row

        # 路径排序与编号：按第一段命中顺序、第二段命中顺序的首次出现次序
        for idx, path in enumerate(paths, start=1):
            path.path_index = idx
            path.fragments.sort(key=lambda f: (f.src_start, f.src_end))
            path.covered2 = sum(f.src_end - f.src_start for f in path.fragments)

        collapsed = len(hits1) - len(reached_hits1)
        collapsed_note = (
            f"；另有 {collapsed} 个第一段候选在第二段无映射，未形成可达路径"
            if collapsed
            else ""
        )

        row.candidates = paths
        if len(paths) > 1:
            row.status = "ambiguous"
            row.reason = (
                f"命中 {len(paths)} 条可达路径（"
                + "、".join(
                    f"路径{p.path_index}: chain {p.chain_id_1}→{p.chain_id_2}"
                    for p in paths
                )
                + "），按要求不自动选择或合并，全部路径列于下方"
                + collapsed_note
            )
            return row

        path = paths[0]
        if path.covered2 == total:
            row.status = "mapped"
            row.reason = ("两段均完全覆盖" + collapsed_note) if collapsed_note else ""
        else:
            row.status = "partial"
            lost1 = total - path.covered1
            lost2 = path.covered1 - path.covered2
            bits = []
            if lost1:
                bits.append(f"第一段缺口/未覆盖区域丢失 {lost1} bp")
            if lost2:
                bits.append(f"第二段缺口/未覆盖区域丢失 {lost2} bp")
            if len(path.fragments) > 1:
                bits.append(f"区间被缺口拆成 {len(path.fragments)} 个有序片段")
            bits.append("部分映射，未映射部分被丢弃，残余片段不视为完全映射")
            row.reason = "；".join(bits) + collapsed_note
        return row

    def convert_all(self, records: list[BedRecord]) -> list[TwoSegRowResult]:
        return [self.convert(r) for r in records]
