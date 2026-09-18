"""Add the Zhou et al. 2026 reference (npj Precision Oncology — 17-method
TME scoring re-analysis) and renumber ALL references by order of first
appearance in the text (Cell Press requirement; the 修回期遗留 item).

Safety design:
- only bracket tokens whose expanded contents are ALL integers within the
  known reference range are treated as citations (CI brackets like
  [0.597–0.947] and text brackets are left untouched);
- ranges are expanded before mapping and re-compressed after mapping;
- post-conditions verified: every new id cited, list count == max id,
  first-appearance order sequential, no marker left unmapped.
A backup of the manuscript is written alongside before any change.
"""
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MS = REPO / "docs" / "manuscript_v34.md"
BAK = REPO / "_process" / "2026-09-17_ops_out" / \
    "manuscript_before_renumber.md"

ZHOU_ENTRY = ("Zhou, Q., Kirshtein, A. and Shahriyari, L. (2026). Towards "
              "the tumor microenvironment scoring methods for immune "
              "checkpoint inhibitor response. *npj Precis. Oncol.*, 10, 88. "
              "doi:10.1038/s41698-025-01221-z")
ZHOU_ID = 99  # temporary id, remapped by the renumbering pass

OLD_MAX = 34


def expand(token):
    """'[3–6]'/'[9,10]'/'[15]' -> [3,4,5,6] / [9,10] / [15]."""
    parts = re.split(r"[,;\u2013-]", token)
    nums = []
    for p in parts:
        p = p.strip()
        if not re.fullmatch(r"\d+", p):
            return None
        nums.append(int(p))
    if len(nums) == 2 and len(token.split("\u2013")) == 2:
        a, b = nums
        if a < b:
            return list(range(a, b + 1))
    if len(nums) >= 2 and len(set(nums)) > 1 and any(
            nums[i + 1] - nums[i] not in (1, -1) for i in range(len(nums) - 1)):
        return None  # ambiguous non-contiguous dash token
    return nums


def compress(nums):
    """[3,4,5,9] -> '3\u20135,9'."""
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:
            out.append(f"{nums[i]}\u2013{nums[j]}")
        elif j == i:
            out.append(str(nums[i]))
        else:
            out.append(f"{nums[i]},{nums[j]}")
        i = j + 1
    return ",".join(out)


def main():
    t = MS.read_text(encoding="utf-8")
    shutil.copy2(MS, BAK)

    i0 = t.find("## REFERENCES")
    head, ref_block = t[:i0], t[i0:]

    # ---- parse the reference list ----
    entries = {}
    for m in re.finditer(r"^(\d+)\.\s+(.+)$", ref_block, re.M):
        entries[int(m.group(1))] = m.group(2).rstrip()
    assert max(entries) == OLD_MAX and len(entries) == OLD_MAX, \
        f"unexpected ref list: {len(entries)} entries, max {max(entries)}"
    entries[ZHOU_ID] = ZHOU_ENTRY

    # ---- insert the in-text Zhou citation ----
    old_intro = ("and EXPRESSO [18] 20 signatures across 3,729 patients "
                 "(Supplementary Table S2 compares their published AUROCs "
                 "with ours)")
    new_intro = ("and EXPRESSO [18] 20 signatures across 3,729 patients; a "
                 "re-analysis of 17 tumour-microenvironment scoring methods "
                 "found none reliable across cancer types [99] "
                 "(Supplementary Table S2 compares their published AUROCs "
                 "with ours)")
    assert old_intro in head, "intro anchor sentence not found"
    head = head.replace(old_intro, new_intro)

    # ---- collect first-appearance order (range-expanded) ----
    order = []
    for m in re.finditer(r"\[([^\[\]]+)\]", head):
        nums = expand(m.group(1))
        if nums is None:
            continue
        if not all(1 <= n <= OLD_MAX or n == ZHOU_ID for n in nums):
            continue
        for n in nums:
            if n not in order:
                order.append(n)
    assert set(order) == set(entries.keys()), \
        f"citation/list mismatch: missing {set(entries) - set(order)}"
    mapping = {old: new for new, old in enumerate(order, 1)}
    print(f"first-appearance order ({len(order)} refs): {order}")

    # ---- rewrite body citations ----
    def rewrite(m):
        nums = expand(m.group(1))
        if nums is None or not all(n in mapping for n in nums):
            return m.group(0)
        return "[" + compress([mapping[n] for n in nums]) + "]"

    head = re.sub(r"\[([^\[\]]+)\]", rewrite, head)

    # ---- rebuild the reference list ----
    lines = ["## REFERENCES", ""]
    for new_id in range(1, len(order) + 1):
        old_id = order[new_id - 1]
        lines.append(f"{new_id}. {entries[old_id]}")
    new_ref_block = "\n".join(lines) + "\n"

    MS.write_text(head + new_ref_block, encoding="utf-8")

    # ---- post-conditions ----
    t2 = MS.read_text(encoding="utf-8")
    body2 = t2.split("## REFERENCES")[0]
    cited = []
    for m in re.finditer(r"\[([^\[\]]+)\]", body2):
        nums = expand(m.group(1))
        if nums and all(1 <= n <= len(order) for n in nums):
            for n in nums:
                if n not in cited:
                    cited.append(n)
    ok_seq = cited == list(range(1, len(order) + 1))
    n_entries = len(re.findall(r"^\d+\. ", t2.split("## REFERENCES")[1],
                               re.M))
    print(f"post: cited sequential={ok_seq}, entries={n_entries}, "
          f"max cited={max(cited)}")
    if not (ok_seq and n_entries == len(order)):
        print("POST-CONDITION FAILED — restoring backup")
        shutil.copy2(BAK, MS)
        return 1
    print("renumber complete: 34 ->", len(order), "references")
    return 0


if __name__ == "__main__":
    sys.exit(main())
