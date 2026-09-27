"""尺寸链约束求解：差分约束系统 / 负权闭环 / 最紧间隙界。

建模（坐标 x_v 表示基准点 v 的坐标，单位微米）：
    约束 m: low_m <= x[to] - x[from] <= high_m
    两条有向边（差分约束标准形式 x[v] - x[u] <= w 即边 u -> v 权 w）：
      from -> to, 权 high_m  （正向：采用上限）
      to   -> from, 权 -low_m（逆向：采用下限）

- 可行性：无负环（Bellman-Ford，N-1 轮松弛 + 前驱链回溯找负权闭环）。
- 最紧界（一致时）：对任意 p, q，x[q]-x[p] 恰为 dist[p][q]（Floyd 全源最短路）。
  故 gap = x[q]-x[p] 的最紧区间为 [ -dist[q][p], dist[p][q] ]；
  任一侧不可达则该点对在约束图中不连通，间隙无界，不判合格。
- 容差预算：每条区间按 [low-k*s, high+k*s] 随预算 k 放宽（s 为该测量步长，
  0 即锁定），边权随 k 单调不减，一致性关于 k 单调，二分求首次一致的最小 k。
- 名义坐标拟合：锚点代表绝对零位（首点坐标固定为 0），|x_i − 名义_i| ≤ D
  化作一对锚点界边；D 的可行性即负环判定、关于 D 单调，二分求最小 D；
  再按录入顺序逐点固定为当前最紧下界 −dist(i, 锚点)，得字典序最小解。
- 垫片组合选配：垫片不改变测量，只把目标点对的有界最紧区间 [L, U] 整体平移
  为 [L+S, U+S]（S 为各安装位“方向 × 所选厚度”之和）。落入合格区间 [T₀, T₁]
  等价于 S ∈ [T₀−L, T₁−U]；按“有符号厚度和 → (最低成本, 最小字典序秩)”
  做分层动态规划（同和状态合并、按后缀可达区间剪枝；每层按预计可达状态数
  在稀疏字典与稠密数组间重选表示，稀疏可达和绝不按整个整数跨度开数组），
  不逐一枚举 ∏选项数 个组合。
- 尺寸链改动影响：同一组基准点（名称与顺序相同）下两版测量清单各自判负环，
  任一版矛盾即只返回各矛盾版本的负环；均一致时对两版分别求 Floyd 全源最短路，
  按基准点录入顺序逐有序异点对比较最紧区间与无界状态，仅报告发生变化的点对，
  并按区间包含关系分为收紧 / 放宽 / 混合，单独区分变为有界与变为无界。
- 矛盾测量最小撤回：只允许删除可撤回测量（至多 16 条），求使剩余清单一致的
  最少删除集合——即原图全部负权闭环的最小碰集；搜索由负环引导分支（只沿当前
  负环上的可撤回编号分支）并以“可撤回编号两两不交的负环打包”作下界剪枝，
  不枚举全部删除子集；并列时按原录入顺序取删除序号列表字典序最小者；
  锁定测量自身矛盾时不可修复，返回锁定部分的负环。
- 测量复核排程：从当前清单选取 1～16 条测量，各设复核耗时、独立不合格概率
  （‰）与先决依赖（只能引用所选测量，成环报错）；复核遇首条不合格即停止。
  求满足依赖且预期总耗时最小的完整顺序，并列时按原录入位置的序号向量取
  字典序最小者；子集动态规划（2^n 状态），不枚举全部排列；概率与期望值
  以精确分数计算与返回。
- 目标间隙精度归因：选取目标点对与 1～10 条待评估测量，未选测量固定为底座。
  对底座加每个待评估子集 S 求目标间隙最紧区间宽度 W(S)，按 Shapley 权重
  |S|!(n−|S|−1)!/n! 把宽度缩减 W(S)−W(S∪{i}) 分摊给各测量，贡献之和恰为
  底座与全量的宽度差；宽度在“终端图”上沿子集递归树增量计算，不枚举排列；
  底座不能给出有界最紧区间时不可归因。
- 测量缺失合格概率：选取目标点对、整数合格区间与 1～12 条可能缺失的测量，
  每条给出独立缺失概率 p_i/1000（未选测量始终有效，0‰/1000‰ 即确定事件）。
  完整清单一致时，每种缺失组合对应“底座（未选测量）＋该组合中仍有效的所选
  测量”的子图，目标间隙有界且最紧区间整体落入合格区间才计为通过；组合权重
  为各条 p_i 或 (1000−p_i) 之积（公分母 1000^n），通过/未通过权重和以既约
  分数返回。各子图最短路在终端图上沿子集递归树增量更新，不为每种组合重跑
  全点最短路；完整清单矛盾时只返回负环。
- 极值共现工况：按顺序接收 1～10 个目标，每项为两个不同基准点构成的有序
  间隙及其最紧下界/上界界别（同一有序点对、同一界别不得重复）。对完整清单
  求各目标的最紧整数界，把“达到该界”作为等式条件（一对反向等权边）叠加到
  原测量图上，求仍能同时成立（无负环）的目标数最多的子集，并列时取输入
  序号列表字典序最小者；返回各目标界值、入选与未入选目标，以及满足原测量
  和入选等式的一组整数坐标——每个原始连通分量以最早录入的基准点为零点
  （分量内整体平移不破坏任何差分约束）。任一目标界无界时报明目标序号、
  不生成工况；原清单矛盾时只返回负环。
所有运算均为整数或精确分数，无浮点误差。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations

MAX_POINTS = 30
MIN_POINTS = 2
MAX_CONSTRAINTS = 80
MIN_CONSTRAINTS = 1
LIMIT_RANGE = 10**9  # 上下限绝对值上限
WINDOW_RANGE = 10**9  # 补测读数 r 的取值范围 ±WINDOW_RANGE


@dataclass(frozen=True)
class Edge:
    u: int
    v: int
    w: int
    cid: str
    # forward=True: 起->止，采用上限；False: 止->起，采用下限
    forward: bool


def validate(payload: dict) -> tuple[list[str], list[dict]]:
    """校验入参，返回 (点名列表, 约束字典列表)；失败抛 ValueError(消息列表)。"""
    errors: list[str] = []

    points = payload.get("points")
    if not isinstance(points, list):
        raise ValueError(["points 必须为字符串数组"])
    norm_points: list[str] = []
    for i, p in enumerate(points):
        name = p.strip() if isinstance(p, str) else ""
        if not name:
            errors.append(f"第 {i + 1} 个基准点：名称为空")
        elif name in norm_points:
            errors.append(f"基准点名称重复：{name}")
        norm_points.append(name)
    if not errors:
        if not (MIN_POINTS <= len(norm_points) <= MAX_POINTS):
            errors.append(f"基准点数量须为 {MIN_POINTS}～{MAX_POINTS} 个（当前 {len(norm_points)}）")

    constraints = payload.get("constraints")
    if not isinstance(constraints, list):
        raise ValueError(["constraints 必须为对象数组"])
    norm_constraints: list[dict] = []
    seen_ids: set[str] = set()
    for i, c in enumerate(constraints):
        prefix = f"约束 {i + 1}"
        if not isinstance(c, dict):
            errors.append(f"{prefix}：格式不是对象")
            continue
        cid = c.get("id")
        cid = cid.strip() if isinstance(cid, str) else ""
        if not cid:
            errors.append(f"{prefix}：编号为空")
        elif cid in seen_ids:
            errors.append(f"{prefix}：编号重复（{cid}）")
        seen_ids.add(cid)

        p_from = c.get("from") if isinstance(c.get("from"), str) else None
        p_to = c.get("to") if isinstance(c.get("to"), str) else None
        p_from = p_from.strip() if p_from else ""
        p_to = p_to.strip() if p_to else ""
        if p_from == "":
            errors.append(f"{prefix}（{cid or '无编号'}）：起点为空")
        elif norm_points and p_from not in norm_points:
            errors.append(f"{prefix}（{cid}）：起点 {p_from} 尚未录入")
        if p_to == "":
            errors.append(f"{prefix}（{cid or '无编号'}）：终点为空")
        elif norm_points and p_to not in norm_points:
            errors.append(f"{prefix}（{cid}）：终点 {p_to} 尚未录入")
        if p_from and p_to and p_from == p_to:
            errors.append(f"{prefix}（{cid}）：起点与终点必须不同")

        low = c.get("low")
        high = c.get("high")
        if not isinstance(low, int) or isinstance(low, bool):
            errors.append(f"{prefix}（{cid}）：下限必须为整数")
            low = None
        if not isinstance(high, int) or isinstance(high, bool):
            errors.append(f"{prefix}（{cid}）：上限必须为整数")
            high = None
        if low is not None and abs(low) > LIMIT_RANGE:
            errors.append(f"{prefix}（{cid}）：下限超出 ±{LIMIT_RANGE}")
        if high is not None and abs(high) > LIMIT_RANGE:
            errors.append(f"{prefix}（{cid}）：上限超出 ±{LIMIT_RANGE}")
        if low is not None and high is not None and low > high:
            errors.append(f"{prefix}（{cid}）：下限 {low} 大于上限 {high}")

        norm_constraints.append({"id": cid, "from": p_from, "to": p_to,
                                 "low": low, "high": high})

    if not errors and not (MIN_CONSTRAINTS <= len(norm_constraints) <= MAX_CONSTRAINTS):
        errors.append(f"约束数量须为 {MIN_CONSTRAINTS}～{MAX_CONSTRAINTS} 条（当前 {len(norm_constraints)}）")

    if errors:
        raise ValueError(errors)
    return norm_points, norm_constraints


def build_edges(points: list[str], constraints: list[dict]) -> list[Edge]:
    idx = {name: i for i, name in enumerate(points)}
    edges: list[Edge] = []
    for c in constraints:
        u = idx[c["from"]]
        v = idx[c["to"]]
        edges.append(Edge(u, v, c["high"], c["id"], True))
        edges.append(Edge(v, u, -c["low"], c["id"], False))
    return edges


def find_negative_cycle(n: int, edges: list[Edge]) -> list[int] | None:
    """Bellman-Ford：返回构成负权闭环的边在 edges 中的下标（沿环顺序），无则 None。"""
    dist = [0] * n                 # 虚拟超源：所有点初值 0
    parent_v = [-1] * n            # 最后一次松弛进入该点的边的起点
    parent_e = [-1] * n            # 最后一次松弛进入该点的边下标

    x = -1
    for _ in range(n):
        x = -1
        for ei, e in enumerate(edges):
            nd = dist[e.u] + e.w
            if nd < dist[e.v]:
                dist[e.v] = nd
                parent_v[e.v] = e.u
                parent_e[e.v] = ei
                x = e.v
        if x == -1:
            return None
    if x == -1:
        return None

    # 第 n 轮仍可松弛的点 x，沿前驱回退 n 步后必在负环上
    y = x
    for _ in range(n):
        y = parent_v[y]
        if y < 0:
            return None

    cycle_edges: list[int] = []
    cur = y
    while True:
        ei = parent_e[cur]
        if ei < 0:
            return None
        cycle_edges.append(ei)
        cur = parent_v[cur]
        if cur == y:
            break
        if len(cycle_edges) > n:
            return None

    cycle_edges.reverse()  # 改为沿环方向（从 y 出发）的顺序
    if sum(edges[i].w for i in cycle_edges) < 0:
        return cycle_edges
    return None


def all_pairs_shortest(n: int, edges: list[Edge]) -> list[list[int | None]]:
    """Floyd 全源最短路（仅在无负环时调用）。不可达为 None。"""
    inf = None
    d: list[list[int | None]] = [[inf] * n for _ in range(n)]
    for i in range(n):
        d[i][i] = 0
    for e in edges:
        old = d[e.u][e.v]
        if old is None or e.w < old:
            d[e.u][e.v] = e.w
    for k in range(n):
        dk = d[k]
        for i in range(n):
            dik = d[i][k]
            if dik is None:
                continue
            di = d[i]
            for j in range(n):
                dkj = dk[j]
                if dkj is None:
                    continue
                cand = dik + dkj
                if di[j] is None or cand < di[j]:
                    di[j] = cand
    return d


def describe_cycle(edges: list[Edge], cycle: list[int]) -> dict:
    segments = []
    total = 0
    for ei in cycle:
        e = edges[ei]
        if e.forward:
            seg = {
                "measurementId": e.cid,
                "direction": "forward",
                "fromIndex": e.u,
                "toIndex": e.v,
                "usedBound": "upper",
                "weight": e.w,
            }
        else:
            seg = {
                "measurementId": e.cid,
                "direction": "reverse",
                "fromIndex": e.u,
                "toIndex": e.v,
                "usedBound": "lower",
                "weight": e.w,
            }
        segments.append(seg)
        total += e.w
    return {"segments": segments, "weightSum": total}


def contradiction_payload(points: list[str], edges: list[Edge], cycle: list[int]) -> dict:
    """负权闭环（边下标列表）-> 对外的矛盾证据结构。"""
    info = describe_cycle(edges, cycle)
    return {
        "segments": [
            {
                "measurementId": s["measurementId"],
                "direction": s["direction"],
                "from": points[s["fromIndex"]],
                "to": points[s["toIndex"]],
                "usedBound": s["usedBound"],
                "weight": s["weight"],
            }
            for s in info["segments"]
        ],
        "weightSum": info["weightSum"],
    }


def analyze(points: list[str], constraints: list[dict]) -> dict:
    edges = build_edges(points, constraints)
    n = len(points)
    cycle = find_negative_cycle(n, edges)
    if cycle is not None:
        return {
            "consistent": False,
            "points": points,
            "contradiction": contradiction_payload(points, edges, cycle),
        }

    dist = all_pairs_shortest(n, edges)
    return {
        "consistent": True,
        "points": points,
        "dist": dist,  # None 表示不可达
    }


def classify(target_low: int, target_high: int, low: int, high: int) -> str:
    """目标区间 [target_low, target_high] 与约束最紧区间 [low, high] 的关系。"""
    if high < target_low or target_high < low:
        return "fail"        # 不相交 -> 超差
    if target_low <= low and high <= target_high:
        return "pass"        # 全落入 -> 合格
    return "review"          # 仅相交 -> 待复核


def query(analysis: dict, p_from: str, p_to: str,
          target_low: int, target_high: int) -> dict:
    points = analysis["points"]
    i = points.index(p_from)
    j = points.index(p_to)
    d_ij = analysis["dist"][i][j]
    d_ji = analysis["dist"][j][i]

    result = {
        "from": p_from,
        "to": p_to,
        "targetLow": target_low,
        "targetHigh": target_high,
    }
    if d_ij is None or d_ji is None:
        result.update({"bounded": False})
        return result

    low = -d_ji
    high = d_ij
    result.update({
        "bounded": True,
        "low": low,
        "high": high,
        "verdict": classify(target_low, target_high, low, high),
    })
    return result


# ============================ 测量清单精简 ============================

def _kept_dist(points: list[str], constraints: list[dict],
               kept: list[bool]) -> list[list[int | None]]:
    """对当前保留约束重算全源最短路（调用方须保证系统一致、无负环）。

    删去一条测量即删去它贡献的两条边；这里只对保留的约束建边，
    复用 all_pairs_shortest 的 Floyd 实现。
    """
    edges = build_edges(points, [c for c, k in zip(constraints, kept) if k])
    return all_pairs_shortest(len(points), edges)


def _first_changed_pair(points: list[str],
                        base: list[list[int | None]],
                        trial: list[list[int | None]]) -> dict:
    """按基准点录入顺序找首个界值/无界状态不同的有序点对，返回其删前删后界。

    删边只会放松最短路：上界 dist[p][q] 只能变大或由有界变无界(None)。
    有序点对 p→q 的最紧区间为 [下界, 上界] = [-dist[q][p], dist[p][q]]，
    逐有序对比较矩阵元素即等价于逐对比较下界、上界与无界状态。
    """
    n = len(points)
    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            if base[p][q] != trial[p][q]:
                upper_before = base[p][q]
                upper_after = trial[p][q]
                return {
                    "from": points[p],
                    "to": points[q],
                    "lowerBefore": None if base[q][p] is None else -base[q][p],
                    "lowerAfter": None if trial[q][p] is None else -trial[q][p],
                    "upperBefore": upper_before,
                    "upperAfter": upper_after,
                    "changedSide": "unbounded" if (upper_before is None) != (upper_after is None)
                    or (base[q][p] is None) != (trial[q][p] is None) else "loosened",
                }
    raise RuntimeError("试删后未找到任何变化点对")


def simplify(points: list[str], constraints: list[dict]) -> dict:
    """测量清单精简（系统一致时）。

    按录入顺序逐条试删：每次从当前保留集合重算全部有序点对的最紧上下界，
    仅当每对的界值与无界状态均与原清单相同时才真正删去；后续判断基于已删结果。
    返回保留/可删编号；对每条保留约束给出“再删它”时按基准点录入顺序找到的
    首个变化点对及其删前、删后界。最终保留集合不可再删。
    """
    n = len(points)
    m = len(constraints)

    # 原清单（全部约束）基准矩阵
    base = _kept_dist(points, constraints, [True] * m)

    kept = [True] * m
    removable: list[str] = []
    for t in range(m):
        trial_kept = kept.copy()
        trial_kept[t] = False
        trial = _kept_dist(points, constraints, trial_kept)
        if trial == base:           # 全部有序点对的界值与无界状态均不变
            kept[t] = False
            removable.append(constraints[t]["id"])

    # 对每条最终保留约束计算“再删它”的见证（基于已删结果）
    current = _kept_dist(points, constraints, kept)
    retained: list[dict] = []
    for t in range(m):
        if not kept[t]:
            continue
        trial_kept = kept.copy()
        trial_kept[t] = False
        trial = _kept_dist(points, constraints, trial_kept)
        witness = _first_changed_pair(points, current, trial)
        retained.append({"id": constraints[t]["id"], "witness": witness})

    return {
        "consistent": True,
        "points": points,
        "simplify": {
            "retainedIds": [constraints[t]["id"] for t in range(m) if kept[t]],
            "removableIds": removable,
            "retained": retained,
        },
    }


# ============================ 测量容差预算 ============================

DEFAULT_STEP = 1
MAX_STEP = 1000


def validate_steps(constraints: list[dict], raw_steps) -> dict[str, int]:
    """校验“编号 -> 步长”映射，返回每条约束的有效步长（未列出的编号取 DEFAULT_STEP）。

    步长须为 0～MAX_STEP 的整数；0 表示锁定（区间不随预算 k 放宽）。
    映射中出现未录入的测量编号（多余编号）即报错。
    """
    if raw_steps is None:
        raw_steps = {}
    if not isinstance(raw_steps, dict):
        raise ValueError(["steps 必须为对象：键为测量编号，值为 0～1000 的整数步长"])

    ids = [c["id"] for c in constraints]
    id_set = set(ids)
    errors: list[str] = []
    given: dict[str, int] = {}
    extras: list[str] = []
    for key, val in raw_steps.items():
        cid = key.strip() if isinstance(key, str) else str(key)
        if cid not in id_set:
            extras.append(cid or "（空编号）")
            continue
        if not isinstance(val, int) or isinstance(val, bool):
            errors.append(f"测量 {cid}：步长必须为整数")
        elif not 0 <= val <= MAX_STEP:
            errors.append(f"测量 {cid}：步长须为 0～{MAX_STEP}（当前 {val}）")
        else:
            given[cid] = val
    if extras:
        errors.append("步长映射包含未录入的测量编号：" + "、".join(extras))
    if errors:
        raise ValueError(errors)
    return {cid: given.get(cid, DEFAULT_STEP) for cid in ids}


def _budget_edges(points: list[str], constraints: list[dict],
                  steps: dict[str, int], k: int) -> list[Edge]:
    """预算 k 下的建边：每条 [low, high] 放宽为 [low − k·步长, high + k·步长]。"""
    idx = {name: i for i, name in enumerate(points)}
    edges: list[Edge] = []
    for c in constraints:
        s = steps[c["id"]]
        u = idx[c["from"]]
        v = idx[c["to"]]
        edges.append(Edge(u, v, c["high"] + k * s, c["id"], True))
        edges.append(Edge(v, u, -(c["low"] - k * s), c["id"], False))
    return edges


def tolerance(points: list[str], constraints: list[dict],
              steps: dict[str, int]) -> dict:
    """测量容差预算：求使整张清单首次一致的最小非负整数 k。

    - k = 0 即原清单：一致则直接返回 0；
    - 锁定测量（步长 0）的区间不随 k 变化，若它们自身已含负权闭环，
      则该矛盾对任意 k 成立，返回不可修复及该负环；
    - 否则边权随 k 单调不减，一致性关于 k 单调，二分最小可行 k；
      返回 k−1 下的负环（再小一步仍矛盾的证据）与 k 下每条测量的
      调整后区间（计算值允许超出原始录入范围 ±LIMIT_RANGE）。
    """
    n = len(points)
    base = {"points": points, "steps": steps}

    edges0 = _budget_edges(points, constraints, steps, 0)
    if find_negative_cycle(n, edges0) is None:
        return {**base, "repairable": True, "k": 0}

    # 锁定测量（步长 0）自身矛盾 -> 不可修复
    locked = [c for c in constraints if steps[c["id"]] == 0]
    if locked:
        locked_edges = build_edges(points, locked)
        cycle = find_negative_cycle(n, locked_edges)
        if cycle is not None:
            return {**base, "repairable": False,
                    "contradiction": contradiction_payload(points, locked_edges, cycle)}

    # 上界：任一简单环至多 n 条边、每边权 >= w_min；含非锁定边的环步长和 >= 1，
    # 故 k >= -n * w_min 时所有环权重非负（k = 0 已矛盾，w_min 必为负）。
    w_min = min(e.w for e in edges0)
    lo, hi = 1, max(1, -n * w_min)
    while lo < hi:
        mid = (lo + hi) // 2
        if find_negative_cycle(n, _budget_edges(points, constraints, steps, mid)) is None:
            hi = mid
        else:
            lo = mid + 1
    k = lo

    # k−1 下的负权闭环（二分不变式保证存在）
    edges_prev = _budget_edges(points, constraints, steps, k - 1)
    cycle_prev = find_negative_cycle(n, edges_prev)
    if cycle_prev is None:
        raise RuntimeError("二分不变式被破坏：k−1 下未找到负权闭环")

    adjusted = [
        {
            "id": c["id"],
            "step": steps[c["id"]],
            "low": c["low"],
            "high": c["high"],
            "adjustedLow": c["low"] - k * steps[c["id"]],
            "adjustedHigh": c["high"] + k * steps[c["id"]],
        }
        for c in constraints
    ]
    return {
        **base,
        "repairable": True,
        "k": k,
        "previousContradiction": contradiction_payload(points, edges_prev, cycle_prev),
        "adjusted": adjusted,
    }


# ============================ 名义坐标拟合 ============================

def validate_nominals(points: list[str], raw_nominals) -> dict[str, int]:
    """校验“点名 -> 名义坐标”映射，返回按基准点录入顺序排列的名义坐标。

    名义坐标须恰好覆盖全部点名（缺点名、含未录入点名均报错）；
    值须为 ±LIMIT_RANGE 内的整数。
    """
    if raw_nominals is None:
        raw_nominals = {}
    if not isinstance(raw_nominals, dict):
        raise ValueError(["nominals 必须为对象：键为基准点名，值为整数名义坐标（微米）"])

    point_set = set(points)
    errors: list[str] = []
    given: dict[str, int] = {}
    seen: set[str] = set()
    extras: list[str] = []
    for key, val in raw_nominals.items():
        name = key.strip() if isinstance(key, str) else str(key)
        if name not in point_set:
            extras.append(name or "（空点名）")
            continue
        seen.add(name)
        if not isinstance(val, int) or isinstance(val, bool):
            errors.append(f"基准点 {name}：名义坐标必须为整数")
        elif abs(val) > LIMIT_RANGE:
            errors.append(f"基准点 {name}：名义坐标超出 ±{LIMIT_RANGE}")
        else:
            given[name] = val
    if extras:
        errors.append("名义坐标映射包含未录入的基准点：" + "、".join(extras))
    missing = [p for p in points if p not in seen]
    if missing:
        errors.append("名义坐标映射缺少基准点：" + "、".join(missing))
    if errors:
        raise ValueError(errors)
    return {p: given[p] for p in points}


def _feasible_potential(n: int, edges: list[Edge]) -> list[int]:
    """Bellman-Ford（虚拟超源，初值全 0）；调用方须保证无负环。返回一组可行势。"""
    dist = [0] * n
    for _ in range(n - 1):
        updated = False
        for e in edges:
            nd = dist[e.u] + e.w
            if nd < dist[e.v]:
                dist[e.v] = nd
                updated = True
        if not updated:
            break
    return dist


def _distances_to(n: int, edges: list[Edge], target: int) -> list[int | None]:
    """无负环图中各点到 target 的最短路（反向图上自 target 的 Bellman-Ford）。"""
    dist: list[int | None] = [None] * n
    dist[target] = 0
    for _ in range(n - 1):
        updated = False
        for e in edges:
            dv = dist[e.v]
            if dv is None:
                continue
            cand = dv + e.w
            du = dist[e.u]
            if du is None or cand < du:
                dist[e.u] = cand
                updated = True
        if not updated:
            break
    return dist


def coordinate_fit(points: list[str], constraints: list[dict],
                   nominals: dict[str, int]) -> dict:
    """名义坐标拟合（系统一致时调用；矛盾检测由调用方先行完成）。

    在满足全部测量的整数坐标方案中（首个基准点实际坐标固定为 0）：
      1. 最小化 D = max |x_i − 名义_i|：|x_i − 名义_i| ≤ D 即一对锚点界边，
         给定 D 的可行性 = 无负环，且关于 D 单调，二分求最小可行 D；
      2. 在 D 最小的方案中按基准点录入顺序取坐标向量字典序最小者：
         差分约束对单变量的投影是区间 [−dist(i,锚点), dist(锚点,i)]，
         固定其中任一整数值都保持可行，故逐点取当前最紧下界即可（贪心）。
    返回 D、各点实际坐标与有符号偏差（实际 − 名义）；计算值允许超出 ±LIMIT_RANGE。
    """
    n = len(points)
    a = n  # 锚点下标：代表绝对零位，x[首点] 经两条 0 权边与锚点相等
    base_edges = build_edges(points, constraints)
    anchor_edges = [Edge(a, 0, 0, "", True), Edge(0, a, 0, "", False)]

    def edges_with_d(limit: int) -> list[Edge]:
        edges = list(base_edges) + list(anchor_edges)
        for i, p in enumerate(points):
            nv = nominals[p]
            edges.append(Edge(a, i, nv + limit, "", True))    # x_i − x_a ≤ 名义_i + D
            edges.append(Edge(i, a, limit - nv, "", False))   # x_a − x_i ≤ D − 名义_i
        return edges

    # D 的上界：任取一组可行势（平移使首点为 0）即为一个可行方案
    pot = _feasible_potential(n + 1, base_edges + anchor_edges)
    hi = max(abs(pot[i] - pot[0] - nominals[p]) for i, p in enumerate(points))
    lo = 0
    while lo < hi:
        mid = (lo + hi) // 2
        if find_negative_cycle(n + 1, edges_with_d(mid)) is None:
            hi = mid
        else:
            lo = mid + 1
    limit = lo

    # 字典序最小化：按录入顺序逐点取当前最紧下界并固定（界边保证各点可达锚点）
    edges = edges_with_d(limit)
    coords = [0] * n
    for i in range(n):
        coords[i] = -_distances_to(n + 1, edges, a)[i]
        edges.append(Edge(a, i, coords[i], "", True))
        edges.append(Edge(i, a, -coords[i], "", False))

    max_dev = max(abs(coords[i] - nominals[p]) for i, p in enumerate(points))
    return {
        "consistent": True,
        "points": points,
        "nominals": {p: nominals[p] for p in points},
        "fit": {
            "maxDeviation": max_dev,
            "coordinates": [
                {
                    "point": p,
                    "nominal": nominals[p],
                    "actual": coords[i],
                    "deviation": coords[i] - nominals[p],
                }
                for i, p in enumerate(points)
            ],
        },
    }


# ============================ 补测读数窗口 ============================

def reading_window(points: list[str], constraints: list[dict],
                   from_point: str, to_point: str,
                   target_from: str, target_to: str,
                   target_low: int, target_high: int) -> dict:
    """补测读数窗口（仅在原清单一致时调用；矛盾检测由调用方先行完成）。

    在基准点 from_point → to_point 之间拟补一条精确测量 [r, r]（整数 μm）。
    加入后：
      - 整张清单须仍一致（无负环）；
      - 目标点对 target_from → target_to 的间隙最紧界须有界，且完全落入
        合格区间 [target_low, target_high]。
    在 r ∈ [−WINDOW_RANGE, WINDOW_RANGE] 的整数范围内，求全部满足条件的 r
    所构成的、不可再扩展（极大）且互不相邻的整数闭区间，按起点升序返回。

    精确补测 [r,r] 只新增两条边：a→b 权 r、b→a 权 −r（a、b 为补测起终点）。
    设原清单（无负环）Floyd 全源最短路为 dist，则任意两点新最紧上界：
        dist'[u][v] = min(dist[u][v],
                          dist[u][a] + r + dist[b][v]  （经 a→b 边）,
                          dist[u][b] − r + dist[a][v]) （经 b→a 边）
    （路径至多使用一次新边：任一路径多次使用两条新边即含环；无负环时去环后
      不更长。各项不可达视为 +∞。）于是：

    1) 一致性：新环必含新边。a→b 边参与的负环对应 dist[b][a] + r < 0，
       即 r < −dist[b][a]（下界 L，dist[b][a] 不可达则无此约束）；
       b→a 边参与的负环对应 dist[a][b] − r < 0，即 r > dist[a][b]
       （上界 U，dist[a][b] 不可达则无此约束）。可行 r ∈ [L, U]。
    2) 有界且合格：对目标上界 dist'[p][q] ≤ target_high、目标下界
       −dist'[q][p] ≥ target_low ⇔ dist'[q][p] ≤ −target_low 各禁用一个 r
       区间；另须两种经新边的路径都存在，否则对应新最紧界仍可能无界
       （原界已不可达时，缺一种经新边路径即永久无界 → 禁用全部 r）。
    所有界均为 r 的分段线性函数，禁用集至多是两个整数闭区间/半直线，与
    [−WINDOW_RANGE, WINDOW_RANGE] 交并后取补集即可，全程 O(n³)，不枚举 r。
    """
    n = len(points)
    idx = {name: i for i, name in enumerate(points)}
    a = idx[from_point]
    b = idx[to_point]
    p = idx[target_from]
    q = idx[target_to]

    dist = all_pairs_shortest(n, build_edges(points, constraints))

    # 一致性给出的 r 范围 [lower, upper]（已含 ±WINDOW_RANGE 裁剪）
    lower = -WINDOW_RANGE
    upper = WINDOW_RANGE
    if dist[b][a] is not None:
        lower = max(lower, -dist[b][a])
    if dist[a][b] is not None:
        upper = min(upper, dist[a][b])

    def forbidden(u: int, v: int, bound: int) -> tuple[int | None, int | None]:
        """使 dist'[u][v] > bound 的整数 r 集合，表示为一个闭区间 [f_lo, f_hi]；
        端为 None 表示该侧无界；空集返回 (1, 0)。

        dist'[u][v] = min(原界, α+r, β−r)，其中
            α = dist[u][a] + dist[b][v]（经 a→b），β = dist[u][b] + dist[a][v]（经 b→a）。
        min > bound 须各项同时 > bound；原界已不可达或 > bound 时只需
            α + r > bound ⇔ r ≥ bound−α+1，β − r > bound ⇔ r ≤ β−bound−1，
        两个半直线取交集即单个闭区间。
        """
        t0 = dist[u][v]
        if t0 is not None and t0 <= bound:
            return 1, 0  # 原界已不超过 bound：无禁用值
        has_plus = dist[u][a] is not None and dist[b][v] is not None
        has_minus = dist[u][b] is not None and dist[a][v] is not None
        if not has_plus and not has_minus:
            return None, None  # 目标界无界且任何 r 下都无法建立连接：全部禁用
        f_lo = None
        f_hi = None
        if has_plus:
            f_lo = bound - (dist[u][a] + dist[b][v]) + 1
        if has_minus:
            f_hi = (dist[u][b] + dist[a][v]) - bound - 1
        if f_lo is not None and f_hi is not None and f_lo > f_hi:
            return 1, 0
        return f_lo, f_hi

    # 目标上界超 target_high、目标下界低于 target_low 各自的禁用区间
    forbidden_intervals = [
        forbidden(p, q, target_high),
        forbidden(q, p, -target_low),
    ]

    # 裁剪到 [lower, upper]（None 端为无界），剔除空集，合并相交/相邻区间
    blocked: list[tuple[int, int]] = []
    for f_lo, f_hi in forbidden_intervals:
        if f_lo is not None and f_hi is not None and f_lo > f_hi:
            continue
        x = lower if f_lo is None else max(f_lo, lower)
        y = upper if f_hi is None else min(f_hi, upper)
        if x <= y:
            blocked.append((x, y))
    blocked.sort()
    merged: list[tuple[int, int]] = []
    for x, y in blocked:
        if merged and x <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], y))
        else:
            merged.append((x, y))

    # 合并后的禁用区间在 [lower, upper] 中的补集，即全部极大、互不相邻窗口
    windows: list[dict] = []
    cur = lower
    for x, y in merged:
        if cur <= x - 1:
            windows.append({"low": cur, "high": x - 1})
        cur = max(cur, y + 1)
    if cur <= upper:
        windows.append({"low": cur, "high": upper})

    return {
        "consistent": True,
        "points": points,
        "from": from_point,
        "to": to_point,
        "targetFrom": target_from,
        "targetTo": target_to,
        "targetLow": target_low,
        "targetHigh": target_high,
        "rangeLow": -WINDOW_RANGE,
        "rangeHigh": WINDOW_RANGE,
        "windows": windows,
    }


# ============================ 垫片组合选配 ============================

MIN_POSITIONS = 1
MAX_POSITIONS = 10
MIN_OPTIONS = 1
MAX_OPTIONS = 6
MAX_THICKNESS = 10**6
MAX_SHIM_COST = 1000


def validate_positions(raw_positions) -> list[dict]:
    """校验垫片安装位定义，返回规范化列表；失败抛 ValueError(消息列表)。

    每位：direction ∈ {+1, −1}；options 1～6 个，每个为
    {thickness: 0～MAX_THICKNESS 整数, cost: 0～MAX_SHIM_COST 整数}。
    """
    if not isinstance(raw_positions, list):
        raise ValueError(["安装位 positions 必须为对象数组"])
    if not (MIN_POSITIONS <= len(raw_positions) <= MAX_POSITIONS):
        raise ValueError([f"安装位数量须为 {MIN_POSITIONS}～{MAX_POSITIONS} 个"
                          f"（当前 {len(raw_positions)}）"])

    errors: list[str] = []
    positions: list[dict] = []
    for i, pos in enumerate(raw_positions):
        prefix = f"安装位 {i + 1}"
        if not isinstance(pos, dict):
            errors.append(f"{prefix}：格式不是对象")
            continue
        direction = pos.get("direction")
        if not isinstance(direction, int) or isinstance(direction, bool) or direction not in (1, -1):
            errors.append(f"{prefix}：方向必须为 +1 或 -1（当前 {direction!r}）")

        raw_options = pos.get("options")
        if not isinstance(raw_options, list):
            errors.append(f"{prefix}：选项 options 必须为数组")
            options = []
        else:
            if not (MIN_OPTIONS <= len(raw_options) <= MAX_OPTIONS):
                errors.append(f"{prefix}：选项数量须为 {MIN_OPTIONS}～{MAX_OPTIONS} 个"
                              f"（当前 {len(raw_options)}）")
            options = []
            for j, opt in enumerate(raw_options):
                oprefix = f"{prefix} 第 {j + 1} 个选项"
                if not isinstance(opt, dict):
                    errors.append(f"{oprefix}：格式不是对象")
                    continue
                thickness = opt.get("thickness")
                cost = opt.get("cost")
                if not isinstance(thickness, int) or isinstance(thickness, bool):
                    errors.append(f"{oprefix}：厚度必须为整数")
                elif not 0 <= thickness <= MAX_THICKNESS:
                    errors.append(f"{oprefix}：厚度须为 0～{MAX_THICKNESS} 微米"
                                  f"（当前 {thickness}）")
                if not isinstance(cost, int) or isinstance(cost, bool):
                    errors.append(f"{oprefix}：成本必须为整数")
                elif not 0 <= cost <= MAX_SHIM_COST:
                    errors.append(f"{oprefix}：成本须为 0～{MAX_SHIM_COST}（当前 {cost}）")
                if isinstance(thickness, int) and not isinstance(thickness, bool) \
                        and isinstance(cost, int) and not isinstance(cost, bool):
                    options.append({"thickness": thickness, "cost": cost})
        positions.append({"direction": direction, "options": options})

    if errors:
        raise ValueError(errors)
    return positions


def _suffix_extrema(positions: list[dict], t: int) -> tuple[int, int]:
    """安装位 t..k-1 可叠加的有符号厚度和的最小/最大值；t == k 时为 (0, 0)。"""
    lo = hi = 0
    for pos in positions[t:]:
        d = pos["direction"]
        vals = [d * o["thickness"] for o in pos["options"]]
        lo += min(vals)
        hi += max(vals)
    return lo, hi


def shim_plan(points: list[str], constraints: list[dict],
              target_from: str, target_to: str,
              target_low: int, target_high: int,
              positions: list[dict]) -> dict:
    """垫片组合选配（仅在原清单一致、目标点对有界时调用）。

    垫片不改动任何测量，只把目标点对 target_from → target_to 的原最紧间隙区间
    [L, U]（有界）整体平移为 [L+S, U+S]，其中每位恰选一个选项、
    S = Σ 方向_i × 所选厚度_i。要求平移后整个区间落入合格区间
    [target_low, target_high]，即 S ∈ [target_low−L, target_high−U]。

    优化目标：先最小化所选选项总成本；并列时按安装位顺序取从 1 开始的
    选项序号向量字典序最小者。算法为按“有符号厚度和”分层的动态规划：
    同一和状态只保留 (最低成本, 最小字典序秩)（打包为一个整数编码），
    并按后缀可达区间剪掉绝无可能落入所需区间的状态；秩是选项序号向量的
    6 进制编码，末层最优状态的秩逐位解码即得选项序号向量。每层按
    “预计可达状态数 / 保留跨度”在稀疏字典与稠密区间数组（偏移 + INF
    哨兵）之间重新选择表示：稀疏可达和只按状态数分配内存，绝不按整个
    整数跨度开数组。状态数是可达有符号和的个数，不是 ∏选项数，
    故不逐一枚举全部组合。
    """
    n = len(points)
    idx = {name: i for i, name in enumerate(points)}
    p = idx[target_from]
    q = idx[target_to]
    dist = all_pairs_shortest(n, build_edges(points, constraints))

    base: dict = {
        "consistent": True,
        "points": points,
        "targetFrom": target_from,
        "targetTo": target_to,
        "targetLow": target_low,
        "targetHigh": target_high,
    }

    d_pq = dist[p][q]
    d_qp = dist[q][p]
    if d_pq is None or d_qp is None:
        # 目标间隙无界：平移一个无界区间没有意义，无法选配
        base["bounded"] = False
        base["feasible"] = False
        base["reason"] = "unbounded"
        base["message"] = (f"{target_from} → {target_to} 的原间隙最紧区间无界，"
                           "垫片只能整体平移区间、不能使其变界，无法选配。")
        return base

    low = -d_qp          # 原最紧下界 L
    high = d_pq          # 原最紧上界 U
    req_lo = target_low - low    # 所需净补偿量下界 A
    req_hi = target_high - high  # 所需净补偿量上界 B
    base.update({
        "bounded": True,
        "originalLow": low,
        "originalHigh": high,
        "requiredLow": req_lo,
        "requiredHigh": req_hi,
    })
    if req_lo > req_hi:
        base["feasible"] = False
        base["reason"] = "targetTooNarrow"
        base["message"] = (
            f"合格区间 [ {target_low} , {target_high} ] μm 容不下原最紧区间"
            f" [ {low} , {high} ] μm（长度 {high - low} μm），"
            "整体平移无法使其落入，垫片选配无解。")
        return base

    k = len(positions)
    # 前向分层 DP，状态键为当前有符号厚度和。状态值把 (最低成本, 最小字典序秩)
    # 打包为一个整数：编码 = 成本·RANK_MOD + 秩（秩 < 6^k = RANK_MOD，无进位），
    # 整数大小序即 (成本, 秩) 的字典序，同和合并与末层取最优都退化为整数比较。
    # 秩把已选选项序号（0 起）按安装位顺序编码为 6 进制数：处理第 t 位时
    # 加上 j·6^(k-1-t)。任意两条前 t 位选择不同的路径，其秩差至少 6^(k-t)，
    # 大于后缀 k−t 位所能贡献的最大增量 6^(k-t)−1，故“同成本取秩小”在
    # 中间层合并状态时与最终的选项序号向量字典序完全一致（后缀选择只依赖
    # 当前和，不依赖历史，故同和的劣势状态可安全丢弃）。剪枝：只保留
    # “后缀可达区间与所需区间有交”的状态（前缀 s 经后缀覆盖
    # [s+suf_lo, s+suf_hi]，与 [req_lo, req_hi] 不交即绝无可能命中）。
    # 层表示逐层重选、可双向切换：
    #   稀疏层：字典 {和: 编码}，只存实际可达的和，内存与状态数成正比；
    #   稠密层：偏移 + 一条平行数组（INF 哨兵表示该和不可达），下标即和，
    #           省去哈希开销，仅当保留跨度不超过 SPAN_LIMIT 且预计可达状态
    #           足够密（不少于跨度的 1/DENSE_FILL）时才使用。
    # 稀疏可达和绝不按整个整数跨度开数组（如选项厚度都是 2×10⁵ 的倍数时，
    # 保留跨度可达 10⁷ 而实际可达和仅数十个）；状态数上界是可达有符号和的
    # 个数而非 ∏选项数，故也不逐一枚举全部组合。
    INF = 1 << 60
    SPAN_LIMIT = 2_000_000
    DENSE_FILL = 8
    RANK_MOD = 6 ** k

    # 前缀位 0..t-1 有符号和的极值（纯函数，按位预算）
    pre_lo = [0] * (k + 1)
    pre_hi = [0] * (k + 1)
    for t, pos in enumerate(positions):
        d = pos["direction"]
        vals = [d * o["thickness"] for o in pos["options"]]
        pre_lo[t + 1] = pre_lo[t] + min(vals)
        pre_hi[t + 1] = pre_hi[t] + max(vals)

    suf0_lo, suf0_hi = _suffix_extrema(positions, 0)
    # 当前层：稀疏 cur_map（{和: 编码}）或稠密 cur_dense（(下标偏移, [编码])）
    # 二者恰有一个非 None；cur_count 为当前状态数（稠密层为非 INF 槽数），
    # 供下一层估计可达状态密度
    cur_map: dict[int, int] | None = (
        {0: 0} if suf0_lo <= req_hi and suf0_hi >= req_lo else {})
    cur_dense: tuple[int, list[int]] | None = None
    cur_count = len(cur_map)
    for t, pos in enumerate(positions):
        d = pos["direction"]
        s_lo, s_hi = _suffix_extrema(positions, t + 1)
        shift = 6 ** (k - t - 1)
        options = pos["options"]
        m = len(options)
        vals = [d * o["thickness"] for o in options]

        # 新层“保留和”区间：候选和理论区间（前缀整体区间 + 每位选项）与
        # 剪枝窗口 [req_lo−s_hi, req_hi−s_lo] 取交；可达和必落在该交集中，
        # 故稀疏、稠密两种表示下保留的状态集合一致
        span_min = pre_lo[t] + min(vals)
        span_max = pre_hi[t] + max(vals)
        keep_lo = max(span_min, req_lo - s_hi)
        keep_hi = min(span_max, req_hi - s_lo)
        if keep_lo > keep_hi:
            cur_map, cur_dense, cur_count = {}, None, 0   # 剪枝后已无可达状态
            continue
        width = keep_hi - keep_lo + 1
        # 逐层重选表示：保留跨度不大、且预计可达状态数（上界为
        # min(跨度, 当前状态数 × 本层选项数)）足够密集时才用稠密数组，
        # 否则用稀疏字典——稀疏大跨度时绝不按整个整数跨度开数组
        use_dense = (width <= SPAN_LIMIT
                     and min(width, cur_count * m) * DENSE_FILL >= width)

        # 上一层的 (前缀和, 编码) 迭代器：稠密源跳过 INF 槽，并按新保留
        # 区间裁剪扫描下标（区间外的和经本层任何选项都到不了保留区间），
        # 避免逐槽遍历
        if cur_dense is not None:
            old_base, ovals = cur_dense
            src_lo = max(0, keep_lo - max(vals) - old_base)
            src_hi = min(len(ovals) - 1, keep_hi - min(vals) - old_base)

            def state_iter(_lo=src_lo, _hi=src_hi):
                for si in range(_lo, _hi + 1):
                    v0 = ovals[si]
                    if v0 != INF:
                        yield old_base + si, v0
        else:
            def state_iter():
                yield from cur_map.items()

        if use_dense:
            offset = keep_lo
            avals = [INF] * width
            count = 0
            for s, v0 in state_iter():
                for j in range(m):
                    ni = s + vals[j] - offset
                    if ni < 0 or ni >= width:
                        continue
                    nv = v0 + options[j]["cost"] * RANK_MOD + j * shift
                    if nv < avals[ni]:
                        if avals[ni] == INF:
                            count += 1
                        avals[ni] = nv
            cur_map, cur_dense, cur_count = None, (offset, avals), count
        else:
            # 稀疏字典：只存实际可达的和，内存与状态数成正比
            nxt: dict[int, int] = {}
            for s, v0 in state_iter():
                for j in range(m):
                    ns = s + vals[j]
                    if ns + s_lo > req_hi or ns + s_hi < req_lo:
                        continue
                    nv = v0 + options[j]["cost"] * RANK_MOD + j * shift
                    old = nxt.get(ns)
                    if old is None or nv < old:
                        nxt[ns] = nv
            cur_map, cur_dense, cur_count = nxt, None, len(nxt)

    best_state = best_val = None
    if cur_dense is not None:
        offset, avals = cur_dense
        i0 = max(0, req_lo - offset)
        i1 = min(len(avals) - 1, req_hi - offset)
        for ni in range(i0, i1 + 1):
            v = avals[ni]
            if v != INF and (best_val is None or v < best_val):
                best_state, best_val = offset + ni, v
    else:
        for s, v in cur_map.items():
            if req_lo <= s <= req_hi and (best_val is None or v < best_val):
                best_state, best_val = s, v
    if best_state is None:
        base["feasible"] = False
        base["reason"] = "noCombination"
        base["message"] = (
            "没有任何一种垫片组合能使净补偿量落入所需区间"
            f" [ {req_lo} , {req_hi} ] μm（原最紧区间整体平移后落入"
            f" [ {target_low} , {target_high} ] μm 所要求的净补偿量），垫片选配无解。")
        return base

    # 编码解包：成本 = 编码 // RANK_MOD，秩 = 编码 % RANK_MOD；
    # 选项序号向量即秩的 6 进制各位（每位 0～选项数−1 < 6）
    best_cost, best_rank = divmod(best_val, RANK_MOD)
    choice = [(best_rank // 6 ** (k - 1 - t)) % 6 for t in range(k)]
    net = best_state
    adjusted_low = low + net
    adjusted_high = high + net

    selections = []
    for t, pos in enumerate(positions):
        j = choice[t]
        opt = pos["options"][j]
        selections.append({
            "position": t + 1,
            "direction": pos["direction"],
            "optionIndex": j + 1,             # 从 1 开始的选项序号
            "thickness": opt["thickness"],
            "cost": opt["cost"],
            "signedThickness": pos["direction"] * opt["thickness"],
            "optionCount": len(pos["options"]),
        })

    base.update({
        "feasible": True,
        "plan": {
            "netCompensation": net,
            "adjustedLow": adjusted_low,
            "adjustedHigh": adjusted_high,
            "totalCost": best_cost,
            "choiceVector": [j + 1 for j in choice],
            "selections": selections,
        },
        "message": (f"已选垫片组合：净补偿量 {net:+d} μm，调整后区间 "
                    f"[ {adjusted_low} , {adjusted_high} ] μm 完全落入合格区间"
                    f" [ {target_low} , {target_high} ] μm，总成本 {best_cost}。"),
    })
    return base


# ============================ 尺寸链改动影响 ============================

def classify_impact(lower_before: int | None, upper_before: int | None,
                    lower_after: int | None, upper_after: int | None) -> str:
    """前后最紧区间（None 表示该侧无界）的变化类别；调用方须保证前后不完全相同。

    有界↔无界的跃迁优先标出：旧有界新无界 -> becameUnbounded（变为无界），
    旧无界新有界 -> becameBounded（变为有界）；其余按区间包含关系（无界侧
    视为 ±∞）：新区间真包含于旧区间 -> tightened（收紧），旧区间真包含于
    新区间 -> loosened（放宽），互不包含 -> mixed（混合）。
    """
    old_bounded = lower_before is not None and upper_before is not None
    new_bounded = lower_after is not None and upper_after is not None
    if old_bounded and not new_bounded:
        return "becameUnbounded"
    if new_bounded and not old_bounded:
        return "becameBounded"

    def le_lower(a: int | None, b: int | None) -> bool:
        """下界 a ≤ b（None 视为 −∞）。"""
        if a is None:
            return True
        if b is None:
            return False
        return a <= b

    def le_upper(a: int | None, b: int | None) -> bool:
        """上界 a ≤ b（None 视为 +∞）。"""
        if b is None:
            return True
        if a is None:
            return False
        return a <= b

    # 新区间 ⊆ 旧区间 ⇔ 旧下界 ≤ 新下界 且 新上界 ≤ 旧上界（前后不同即真包含）
    if le_lower(lower_before, lower_after) and le_upper(upper_after, upper_before):
        return "tightened"
    if le_lower(lower_after, lower_before) and le_upper(upper_before, upper_after):
        return "loosened"
    return "mixed"


def chain_impact(points: list[str], base_constraints: list[dict],
                 draft_constraints: list[dict]) -> dict:
    """尺寸链改动影响：同一组基准点（名称与顺序相同）下基线与编辑稿的对比。

    两版各自判负权闭环：任一版矛盾即返回各矛盾版本的负环
    （baselineContradiction / draftContradiction），不生成点对差异。
    两版均一致时，按基准点录入顺序逐有序异点对 (p→q, p≠q) 比较最紧区间
    [−dist[q][p], dist[p][q]] 与无界状态，仅返回发生变化的点对（changes，
    含前后界——无界为 None——与变化类别 kind）。
    """
    n = len(points)
    base_edges = build_edges(points, base_constraints)
    draft_edges = build_edges(points, draft_constraints)
    base_cycle = find_negative_cycle(n, base_edges)
    draft_cycle = find_negative_cycle(n, draft_edges)
    if base_cycle is not None or draft_cycle is not None:
        result: dict = {"consistent": False, "points": points}
        if base_cycle is not None:
            result["baselineContradiction"] = \
                contradiction_payload(points, base_edges, base_cycle)
        if draft_cycle is not None:
            result["draftContradiction"] = \
                contradiction_payload(points, draft_edges, draft_cycle)
        return result

    base_dist = all_pairs_shortest(n, base_edges)
    draft_dist = all_pairs_shortest(n, draft_edges)
    changes: list[dict] = []
    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            lower_before = None if base_dist[q][p] is None else -base_dist[q][p]
            upper_before = base_dist[p][q]
            lower_after = None if draft_dist[q][p] is None else -draft_dist[q][p]
            upper_after = draft_dist[p][q]
            if lower_before == lower_after and upper_before == upper_after:
                continue
            changes.append({
                "from": points[p],
                "to": points[q],
                "lowerBefore": lower_before,
                "upperBefore": upper_before,
                "lowerAfter": lower_after,
                "upperAfter": upper_after,
                "kind": classify_impact(lower_before, upper_before,
                                        lower_after, upper_after),
            })
    return {"consistent": True, "points": points, "changes": changes}


# ============================ 矛盾测量最小撤回 ============================

MAX_REMOVABLE = 16


def validate_removable(constraints: list[dict], raw_ids) -> list[str]:
    """校验可撤回编号列表，返回按录入顺序排列的编号；失败抛 ValueError(消息列表)。

    列表缺省（None）按空列表处理；元素须为已录入的测量编号（未知编号报错）、
    不得重复，去重后至多 MAX_REMOVABLE 条。
    """
    if raw_ids is None:
        return []
    if not isinstance(raw_ids, list):
        raise ValueError(["removableIds 必须为数组：元素为已录入的测量编号"])

    id_set = {c["id"] for c in constraints}
    errors: list[str] = []
    seen: set[str] = set()
    unknown: list[str] = []
    for j, rid in enumerate(raw_ids):
        cid = rid.strip() if isinstance(rid, str) else ""
        if not cid:
            errors.append(f"可撤回列表第 {j + 1} 项：编号为空或不是字符串")
            continue
        if cid not in id_set:
            unknown.append(cid)
            continue
        if cid in seen:
            errors.append(f"可撤回编号重复：{cid}")
            continue
        seen.add(cid)
    if unknown:
        errors.append("可撤回列表包含未录入的测量编号：" + "、".join(unknown))
    if len(seen) > MAX_REMOVABLE:
        errors.append(f"可撤回测量至多 {MAX_REMOVABLE} 条（当前 {len(seen)} 条）")
    if errors:
        raise ValueError(errors)
    order = {c["id"]: i for i, c in enumerate(constraints)}
    return sorted(seen, key=lambda cid: order[cid])


def min_removal(points: list[str], constraints: list[dict],
                removable_ids: list[str]) -> dict:
    """矛盾测量最小撤回：只允许删除可撤回测量，求使剩余清单一致的最少删除。

    - 原清单一致：返回空撤回集（removedIds 为空）；
    - 锁定测量（未标记可撤回）自身已含负权闭环：删除可撤回测量无法消除
      该闭环，返回不可修复及仅由锁定测量构成的负环；
    - 否则求使剩余清单一致的最少删除集合；并列时按原录入顺序取删除序号
      （从 1 开始的录入位置）列表字典序最小者。

    删除集合 D 可行 ⇔ 原图的每条负权闭环都至少含 D 中一条测量（最小碰集）。
    搜索不枚举删除子集：每次在当前图中找一条负环，只需删除环上某一条可撤回
    测量即可打破该环，故只沿环上的可撤回编号分支；并以“可撤回编号两两不交
    的负环打包”作为下界剪枝（打包出的每个环都需要一次互不相同的删除）。
    删除全部可撤回测量后只剩锁定测量（已验证一致），故解必然存在。
    """
    n = len(points)
    m = len(constraints)
    index_of = {c["id"]: i for i, c in enumerate(constraints)}
    removable_idx = sorted(index_of[cid] for cid in removable_ids)
    removable_set = set(removable_idx)
    base = {"points": points,
            "removableIds": [constraints[i]["id"] for i in removable_idx]}

    edges_full = build_edges(points, constraints)
    cycle_full = find_negative_cycle(n, edges_full)
    if cycle_full is None:
        # 原清单一致：空撤回集
        return {**base, "consistent": True, "removal": {
            "removedIds": [], "removedIndices": [], "removedCount": 0,
            "remainingIds": [c["id"] for c in constraints],
        }}
    original = contradiction_payload(points, edges_full, cycle_full)

    locked = [c for i, c in enumerate(constraints) if i not in removable_set]
    locked_edges = build_edges(points, locked)
    locked_cycle = find_negative_cycle(n, locked_edges)
    if locked_cycle is not None:
        # 锁定测量自身矛盾：任何撤回都无法修复
        return {**base, "consistent": False, "repairable": False,
                "originalContradiction": original,
                "contradiction": contradiction_payload(points, locked_edges,
                                                       locked_cycle)}

    # build_edges 按约束顺序每条追加两条边：约束 i 的边即 edges_full[2i:2i+2]
    cons_edges = [edges_full[2 * i: 2 * i + 2] for i in range(m)]

    def current_edges(deleted: frozenset) -> list[Edge]:
        es: list[Edge] = []
        for i in range(m):
            if i not in deleted:
                es += cons_edges[i]
        return es

    def packing_lower_bound(deleted: frozenset) -> int:
        """还需删除条数的下界：可撤回编号两两不交的负环个数。

        在“删去 deleted、且候选池外的可撤回约束也视为删去”的图中反复找负环，
        每找到一个就把环上可撤回约束移出候选池；这些环的可撤回集合两两不交，
        每个环都至少需要一次互不相同的删除。锁定部分一致（上方已验证），
        故找到的环上必有可撤回约束，循环必终止。
        """
        lb = 0
        pool = removable_set - deleted
        while True:
            es: list[Edge] = []
            for i in range(m):
                if i in deleted or (i in removable_set and i not in pool):
                    continue
                es += cons_edges[i]
            cyc = find_negative_cycle(n, es)
            if cyc is None:
                return lb
            lb += 1
            pool -= {index_of[es[ei].cid] for ei in cyc}

    best: list[int] = []

    def dfs(deleted: frozenset) -> None:
        nonlocal best
        if best and len(deleted) + packing_lower_bound(deleted) > len(best):
            return  # 下界已超过已知最优解
        es = current_edges(deleted)
        cyc = find_negative_cycle(n, es)
        if cyc is None:
            cand = sorted(deleted)
            if not best or (len(cand), cand) < (len(best), best):
                best = cand
            return
        if best and len(deleted) + 1 > len(best):
            return  # 再删一条即超过已知最优解
        # 只沿该负环上的可撤回编号分支（锁定部分一致，环上必有可撤回约束）
        branch = sorted({index_of[es[ei].cid] for ei in cyc}
                        & (removable_set - deleted))
        for i in branch:
            dfs(deleted | {i})

    dfs(frozenset())
    removed = best  # 删除全部可撤回测量必可行，best 必已赋值
    removed_set = set(removed)
    return {**base, "consistent": False, "repairable": True,
            "originalContradiction": original,
            "removal": {
                "removedIds": [constraints[i]["id"] for i in removed],
                "removedIndices": [i + 1 for i in removed],
                "removedCount": len(removed),
                "remainingIds": [c["id"] for i, c in enumerate(constraints)
                                 if i not in removed_set],
            }}


# ============================ 测量复核排程 ============================

MIN_INSPECTION_ITEMS = 1
MAX_INSPECTION_ITEMS = 16
MIN_REVIEW_DURATION = 1
MAX_REVIEW_DURATION = 1000
MIN_FAIL_PERMILLE = 0
MAX_FAIL_PERMILLE = 999  # 不合格概率（‰）；合格概率 = (1000 − p) / 1000 ≥ 1/1000


def _find_dep_cycle(items: list[dict]) -> list[str] | None:
    """在“必须先复核”依赖关系上找一条环（用于报错）；无环返回 None。

    边 X → D 表示 X 依赖 D（D 必须先复核）；有向环即互相等待的依赖循环。
    """
    deps_of = {it["id"]: it["deps"] for it in items}
    color: dict[str, int] = {}   # 缺省 0 未访 / 1 在栈上 / 2 已完成
    stack: list[str] = []
    found: list[str] = []

    def dfs(u: str) -> bool:
        color[u] = 1
        stack.append(u)
        for v in deps_of[u]:
            if color.get(v, 0) == 0:
                if dfs(v):
                    return True
            elif color[v] == 1:
                found.extend(stack[stack.index(v):])
                return True
        stack.pop()
        color[u] = 2
        return False

    for it in items:
        if color.get(it["id"], 0) == 0 and dfs(it["id"]):
            return found
    return None


def validate_inspection_items(constraints: list[dict], raw_items) -> list[dict]:
    """校验复核排程的测量选择与排程参数，返回按录入位置排序的条目列表。

    每条：{"id", "entryIndex"(从 1 起的录入位置), "duration", "failPermille", "deps"}。
    规则：从当前清单选取 1～16 条测量（须已录入、不得重复）；复核耗时为
    1～1000 的整数；独立不合格概率为 0～999 的整数（‰）；依赖只能引用所选
    测量（不得引用未选取的测量、不得重复、不得依赖自身）；依赖成环时报错。
    """
    if not isinstance(raw_items, list):
        raise ValueError(["items 必须为数组：元素为 {id, duration, failPermille, deps} 对象"])
    if not (MIN_INSPECTION_ITEMS <= len(raw_items) <= MAX_INSPECTION_ITEMS):
        raise ValueError([f"复核测量须选取 {MIN_INSPECTION_ITEMS}～{MAX_INSPECTION_ITEMS} 条"
                          f"（当前 {len(raw_items)} 条）"])

    entry_of = {c["id"]: i + 1 for i, c in enumerate(constraints)}  # 录入位置（1 起）
    errors: list[str] = []
    items: list[dict] = []
    seen: set[str] = set()
    for k, raw in enumerate(raw_items):
        prefix = f"复核条目 {k + 1}"
        if not isinstance(raw, dict):
            errors.append(f"{prefix}：格式不是对象")
            continue
        cid = raw.get("id")
        cid = cid.strip() if isinstance(cid, str) else ""
        if not cid:
            errors.append(f"{prefix}：编号为空")
        elif cid not in entry_of:
            errors.append(f"{prefix}：测量 {cid} 尚未录入")
        elif cid in seen:
            errors.append(f"{prefix}：测量 {cid} 重复选取")
        else:
            seen.add(cid)

        duration = raw.get("duration")
        if not isinstance(duration, int) or isinstance(duration, bool):
            errors.append(f"{prefix}（{cid or '无编号'}）：复核耗时必须为整数")
            duration = None
        elif not MIN_REVIEW_DURATION <= duration <= MAX_REVIEW_DURATION:
            errors.append(f"{prefix}（{cid}）：复核耗时须为 {MIN_REVIEW_DURATION}～"
                          f"{MAX_REVIEW_DURATION}（当前 {duration}）")
            duration = None

        fail = raw.get("failPermille")
        if not isinstance(fail, int) or isinstance(fail, bool):
            errors.append(f"{prefix}（{cid or '无编号'}）：不合格概率必须为整数（‰）")
            fail = None
        elif not MIN_FAIL_PERMILLE <= fail <= MAX_FAIL_PERMILLE:
            errors.append(f"{prefix}（{cid}）：不合格概率须为 {MIN_FAIL_PERMILLE}～"
                          f"{MAX_FAIL_PERMILLE}‰（当前 {fail}）")
            fail = None

        deps: list[str] = []
        raw_deps = raw.get("deps")
        if raw_deps is None:
            raw_deps = []
        if not isinstance(raw_deps, list):
            errors.append(f"{prefix}（{cid or '无编号'}）：依赖 deps 必须为编号数组")
        else:
            dep_seen: set[str] = set()
            for d in raw_deps:
                did = d.strip() if isinstance(d, str) else ""
                if not did:
                    errors.append(f"{prefix}（{cid or '无编号'}）：依赖编号为空或不是字符串")
                    continue
                if did == cid:
                    errors.append(f"{prefix}（{cid}）：不能依赖自身")
                    continue
                if did in dep_seen:
                    errors.append(f"{prefix}（{cid}）：依赖 {did} 重复")
                    continue
                dep_seen.add(did)
                deps.append(did)
        items.append({"id": cid, "duration": duration, "failPermille": fail, "deps": deps})
    if errors:
        raise ValueError(errors)

    # 依赖只能引用所选测量
    selected = {it["id"] for it in items}
    for k, it in enumerate(items):
        outside = [d for d in it["deps"] if d not in selected]
        if outside:
            errors.append(f"复核条目 {k + 1}（{it['id']}）：依赖只能引用所选测量，"
                          f"未选取：" + "、".join(outside))
    if errors:
        raise ValueError(errors)

    cycle = _find_dep_cycle(items)
    if cycle is not None:
        raise ValueError(["复核依赖成环：" + " → ".join(cycle + [cycle[0]])])

    for it in items:
        it["entryIndex"] = entry_of[it["id"]]
    items.sort(key=lambda it: it["entryIndex"])
    return items


def _fraction_payload(fr: Fraction) -> dict:
    """精确分数 -> 对外结构：既约分子/分母（字符串，可超 JS 安全整数）与精确十进制。

    本模块的分母都是 1000 的幂约分而来（只含因子 2 与 5），十进制展开必有限。
    """
    num, den = fr.numerator, fr.denominator
    d = den
    twos = 0
    while d % 2 == 0:
        d //= 2
        twos += 1
    fives = 0
    while d % 5 == 0:
        d //= 5
        fives += 1
    if d != 1:
        raise RuntimeError("分数分母含 2/5 以外的因子，十进制展开不有限")
    places = max(twos, fives)
    scaled = num * 10 ** places // den
    digits = str(scaled).rjust(places + 1, "0")
    if places:
        ip, fp = digits[:-places], digits[-places:].rstrip("0")
        decimal = ip if not fp else f"{ip}.{fp}"
    else:
        decimal = digits
    return {"num": str(num), "den": str(den), "decimal": decimal}


def inspection_order(points: list[str], constraints: list[dict],
                     items: list[dict]) -> dict:
    """测量复核排程：求满足依赖且预期总耗时最小的完整顺序（子集 DP，不枚举排列）。

    过程：按某一完整顺序逐条复核，遇到首条不合格即停止。顺序 π 的预期总耗时
        E(π) = Σ_k t_{π(k)} · ∏_{j<k} ( s_{π(j)} / 1000 )
    其中 t 为复核耗时、s = 1000 − 不合格概率‰（合格概率的千分分子，s ≥ 1）。
    优化目标：E(π) 最小；并列时按原录入位置的序号向量取字典序最小者。

    子集 DP：令 H(R) = 1000^|R| ·（在“此前各项均合格”条件下处理剩余集合 R 的
    最小预期耗时），则 H(∅) = 0，且
        H(R) = min_{i ∈ R, deps(i) ∩ R = ∅}  t_i·1000^|R| + s_i·H(R \\ {i})
    （i 可作为 R 中下一条 ⇔ 其全部先决测量已不在 R 中；依赖无环保证每个非空
    子集都有可选项）。同一 R 的各候选分母同为 1000^|R|，整数比较即精确比较；
    同值候选取录入位置最小者（首元素即决定序号向量字典序），递归还原即得全局
    字典序最小顺序。状态数 2^n（n ≤ 16）、每状态至多 n 个转移，不枚举 n! 个
    排列；全程整数运算，概率与期望值以既约分数返回。
    """
    n = len(items)
    durations = [it["duration"] for it in items]
    passes = [1000 - it["failPermille"] for it in items]
    entry = [it["entryIndex"] for it in items]
    index_of = {it["id"]: k for k, it in enumerate(items)}
    deps_mask = [0] * n
    for k, it in enumerate(items):
        mask = 0
        for d in it["deps"]:
            mask |= 1 << index_of[d]
        deps_mask[k] = mask

    pow1000 = [1]
    for _ in range(n):
        pow1000.append(pow1000[-1] * 1000)

    size = 1 << n
    H = [0] * size
    choice = [-1] * size
    for mask in range(1, size):
        weight = pow1000[mask.bit_count()]
        best = None
        best_i = -1
        bits = mask
        while bits:
            bit = bits & (-bits)
            bits ^= bit
            i = bit.bit_length() - 1
            if deps_mask[i] & mask:
                continue
            cand = durations[i] * weight + passes[i] * H[mask ^ bit]
            if best is None or cand < best or (cand == best and entry[i] < entry[best_i]):
                best = cand
                best_i = i
        if best_i < 0:
            raise RuntimeError("依赖无环但子集没有可执行的测量")
        H[mask] = best
        choice[mask] = best_i

    order_idx: list[int] = []
    mask = size - 1
    while mask:
        i = choice[mask]
        order_idx.append(i)
        mask ^= 1 << i

    # 每步执行概率（此前各项均合格的概率）与预期总耗时，均为精确分数
    steps: list[dict] = []
    prob = Fraction(1, 1)
    for seq, i in enumerate(order_idx, 1):
        it = items[i]
        steps.append({
            "seq": seq,
            "id": it["id"],
            "entryIndex": it["entryIndex"],
            "duration": it["duration"],
            "failPermille": it["failPermille"],
            "deps": it["deps"],
            "probability": _fraction_payload(prob),
        })
        prob *= Fraction(passes[i], 1000)
    expected = Fraction(H[size - 1], pow1000[n])

    return {
        "points": points,
        "itemCount": n,
        "plan": {
            "order": steps,
            "indexVector": [items[i]["entryIndex"] for i in order_idx],
            "expected": _fraction_payload(expected),
        },
    }


# ============================ 目标间隙精度归因 ============================

MIN_EVALUATE = 1
MAX_EVALUATE = 10


def validate_evaluate_ids(constraints: list[dict], raw_ids) -> list[str]:
    """校验待评估测量编号列表，返回按录入顺序排列的编号；失败抛 ValueError(消息列表)。

    须选取 1～10 条已录入的测量（未知编号报错）、不得重复。
    """
    if not isinstance(raw_ids, list):
        raise ValueError(["measurementIds 必须为数组：元素为待评估的测量编号（1～10 条）"])
    if not (MIN_EVALUATE <= len(raw_ids) <= MAX_EVALUATE):
        raise ValueError([f"待评估测量须选取 {MIN_EVALUATE}～{MAX_EVALUATE} 条"
                          f"（当前 {len(raw_ids)} 条）"])

    id_set = {c["id"] for c in constraints}
    errors: list[str] = []
    seen: set[str] = set()
    unknown: list[str] = []
    for j, rid in enumerate(raw_ids):
        cid = rid.strip() if isinstance(rid, str) else ""
        if not cid:
            errors.append(f"待评估列表第 {j + 1} 项：编号为空或不是字符串")
            continue
        if cid not in id_set:
            unknown.append(cid)
            continue
        if cid in seen:
            errors.append(f"待评估编号重复：{cid}")
            continue
        seen.add(cid)
    if unknown:
        errors.append("待评估列表包含未录入的测量编号：" + "、".join(unknown))
    if errors:
        raise ValueError(errors)
    order = {c["id"]: i for i, c in enumerate(constraints)}
    return sorted(seen, key=lambda cid: order[cid])


def _relax_terminal_edge(d: list[list[int | None]], u: int, v: int, w: int) -> None:
    """向终端图加入一条有向边 u→v（权 w）后就地更新最短路矩阵（None 表示不可达）。

    无负环时任一最短路至多经过新边一次（经过两次即含环，去环后不更长），故
        d'[i][j] = min(d[i][j], d[i][u] + w + d[v][j])。
    """
    t = len(d)
    col_u = [d[i][u] for i in range(t)]
    dv = d[v]
    for i in range(t):
        diu = col_u[i]
        if diu is None:
            continue
        row = d[i]
        via = diu + w
        for j in range(t):
            dvj = dv[j]
            if dvj is None:
                continue
            cand = via + dvj
            cur = row[j]
            if cur is None or cand < cur:
                row[j] = cand


def _subset_widths(base_mat: list[list[int | None]], p_loc: int, q_loc: int,
                   extra: list[list[tuple[int, int, int]]], k: int) -> list[int]:
    """全部 2^k 个待评估子集的目标间隙最紧区间宽度（沿子集递归树增量计算）。

    base_mat 为终端图初始最短路矩阵（底座全源最短路在终端间的汇总，None 不可达）；
    extra[b] 为第 b 条待评估测量贡献的两条有向边（终端下标）。递归树每个节点
    对应一个子集：从父节点矩阵复制后加入一条测量的两条边（各 O(|T|²)），
    整棵树共 2^k 个节点，不枚举 k! 个排列。调用方须保证底座已使目标点对
    连通（故所有超集均有界，宽度必为整数）。
    """
    size = 1 << k
    widths = [0] * size

    def dfs(mask: int, b0: int, d: list[list[int | None]]) -> None:
        widths[mask] = d[p_loc][q_loc] + d[q_loc][p_loc]
        for b in range(b0, k):
            d2 = [row[:] for row in d]
            for u, v, w in extra[b]:
                _relax_terminal_edge(d2, u, v, w)
            dfs(mask | (1 << b), b + 1, d2)

    dfs(0, 0, [row[:] for row in base_mat])
    return widths


def precision_attribution(points: list[str], constraints: list[dict],
                          target_from: str, target_to: str,
                          eval_ids: list[str]) -> dict:
    """目标间隙精度归因（仅在原清单一致时调用；矛盾检测由调用方先行完成）。

    未选取的测量固定为底座；待评估测量共 k 条（1～10）。对底座加每个子集 S
    求目标间隙最紧区间宽度 W(S)（底座即有界，故所有超集均有界），按 Shapley
    权重 |S|!(k−|S|−1)!/k! 把宽度缩减 W(S)−W(S∪{i}) 分摊给测量 i：
        φ_i = Σ_{S⊆N\\{i}}  |S|!(k−|S|−1)!/k! · ( W(S) − W(S∪{i}) )
    贡献之和恰为 W(∅)−W(N)（底座宽度 − 全量宽度，Shapley 效率性）。

    各 W(S) 在“终端图”上增量计算：终端 = 目标点对 + 各待评估测量端点
    （至多 2+2k 个），底座 Floyd 全源最短路汇总为终端间初始边（真实路径按
    终端分段后不短于对应底座最短路，反向展开亦可达，故终端图距离与原图
    一致）；随后沿子集递归树逐条加入待评估测量，每条化作两条有向边按
    _relax_terminal_edge 更新。贡献以既约分数返回（分母可含 3、7 等因子，
    不一定有有限十进制展开），排名按贡献降序、并列按原录入顺序。
    """
    n = len(points)
    idx = {name: i for i, name in enumerate(points)}
    p = idx[target_from]
    q = idx[target_to]
    eval_set = set(eval_ids)
    base_constraints = [c for c in constraints if c["id"] not in eval_set]
    eval_constraints = [c for c in constraints if c["id"] in eval_set]  # 录入顺序

    base: dict = {
        "consistent": True,
        "points": points,
        "targetFrom": target_from,
        "targetTo": target_to,
        "selectedIds": [c["id"] for c in eval_constraints],
    }

    # 底座（未选取测量）全源最短路
    dist0 = all_pairs_shortest(n, build_edges(points, base_constraints))
    if dist0[p][q] is None or dist0[q][p] is None:
        # 底座不能给出有界最紧区间：宽度无从定义，不可归因
        base["attributable"] = False
        base["reason"] = "baseUnbounded"
        base["message"] = (
            f"未选取的测量构成底座后，{target_from} → {target_to} 在约束图中不连通，"
            "最紧间隙区间无界、宽度无从定义，不可归因；请调整待评估选取或目标点对。")
        return base

    # 终端：目标点对 + 各待评估测量端点（去重，保持出现顺序）
    terminals = [p, q]
    seen_t = {p, q}
    for c in eval_constraints:
        for name in (c["from"], c["to"]):
            t = idx[name]
            if t not in seen_t:
                seen_t.add(t)
                terminals.append(t)
    loc = {t: i for i, t in enumerate(terminals)}

    base_mat = [[dist0[u][v] for v in terminals] for u in terminals]
    extra = []
    for c in eval_constraints:
        u = loc[idx[c["from"]]]
        v = loc[idx[c["to"]]]
        extra.append([(u, v, c["high"]), (v, u, -c["low"])])

    k = len(eval_constraints)
    widths = _subset_widths(base_mat, 0, 1, extra, k)

    # Shapley 分摊：φ_i = Σ_{S∌i} |S|!(k−|S|−1)!/k! · (W(S) − W(S∪{i}))
    fact = [1]
    for i in range(1, k + 1):
        fact.append(fact[-1] * i)
    size = 1 << k
    contribs: list[Fraction] = []
    for i in range(k):
        bit = 1 << i
        num = 0
        for mask in range(size):
            if mask & bit:
                continue
            s = mask.bit_count()
            num += fact[s] * fact[k - s - 1] * (widths[mask] - widths[mask | bit])
        contribs.append(Fraction(num, fact[k]))

    width_reduction = widths[0] - widths[size - 1]
    if sum(contribs) != width_reduction:
        raise RuntimeError("Shapley 分摊不变式被破坏：贡献之和不等于底座与全量的宽度差")

    # 排名：贡献降序，并列按原录入顺序（从 1 开始的录入位置）
    entry_of = {c["id"]: i + 1 for i, c in enumerate(constraints)}
    entry = [entry_of[c["id"]] for c in eval_constraints]
    order = sorted(range(k), key=lambda i: (-contribs[i], entry[i]))
    contributions = []
    for r, i in enumerate(order, 1):
        fr = contribs[i]
        contributions.append({
            "id": eval_constraints[i]["id"],
            "entryIndex": entry[i],
            "rank": r,
            "contribution": {"num": str(fr.numerator), "den": str(fr.denominator)},
        })

    base.update({
        "attributable": True,
        "baseWidth": widths[0],
        "fullWidth": widths[size - 1],
        "widthReduction": width_reduction,
        "contributions": contributions,
    })
    return base


# ============================ 测量缺失合格概率 ============================

MIN_RELIABILITY = 1
MAX_RELIABILITY = 12
MIN_MISS_PERMILLE = 0
MAX_MISS_PERMILLE = 1000  # 缺失概率（‰）；0‰=必然有效、1000‰=必然缺失，按确定事件处理


def validate_reliability_items(constraints: list[dict], raw_items) -> list[dict]:
    """校验“可能缺失测量”条目列表，返回按录入顺序排列的条目；失败抛 ValueError(消息列表)。

    每条：{"id", "entryIndex"(从 1 起的录入位置), "missingPermille"}。
    规则：选取 1～12 条已录入测量（未知编号报错）、不得重复；独立缺失概率为
    0～1000 的整数（‰；0‰ 必然有效、1000‰ 必然缺失）。
    """
    if not isinstance(raw_items, list):
        raise ValueError(["items 必须为数组：元素为 {id, missingPermille} 对象"])
    if not (MIN_RELIABILITY <= len(raw_items) <= MAX_RELIABILITY):
        raise ValueError([f"可能缺失的测量须选取 {MIN_RELIABILITY}～{MAX_RELIABILITY} 条"
                          f"（当前 {len(raw_items)} 条）"])

    entry_of = {c["id"]: i + 1 for i, c in enumerate(constraints)}  # 录入位置（1 起）
    errors: list[str] = []
    items: list[dict] = []
    seen: set[str] = set()
    for k, raw in enumerate(raw_items):
        prefix = f"缺失概率条目 {k + 1}"
        if not isinstance(raw, dict):
            errors.append(f"{prefix}：格式不是对象")
            continue
        cid = raw.get("id")
        cid = cid.strip() if isinstance(cid, str) else ""
        if not cid:
            errors.append(f"{prefix}：编号为空")
        elif cid not in entry_of:
            errors.append(f"{prefix}：测量 {cid} 尚未录入")
        elif cid in seen:
            errors.append(f"{prefix}：测量 {cid} 重复选取")
        else:
            seen.add(cid)

        miss = raw.get("missingPermille")
        if not isinstance(miss, int) or isinstance(miss, bool):
            errors.append(f"{prefix}（{cid or '无编号'}）：缺失概率必须为整数（‰）")
            miss = None
        elif not MIN_MISS_PERMILLE <= miss <= MAX_MISS_PERMILLE:
            errors.append(f"{prefix}（{cid}）：缺失概率须为 {MIN_MISS_PERMILLE}～"
                          f"{MAX_MISS_PERMILLE}‰（当前 {miss}）")
            miss = None
        items.append({"id": cid, "missingPermille": miss})
    if errors:
        raise ValueError(errors)

    for it in items:
        it["entryIndex"] = entry_of[it["id"]]
    items.sort(key=lambda it: it["entryIndex"])
    return items


def qualification_reliability(points: list[str], constraints: list[dict],
                              target_from: str, target_to: str,
                              target_low: int, target_high: int,
                              items: list[dict]) -> dict:
    """测量缺失合格概率（仅在完整清单一致时调用；矛盾检测由调用方先行完成）。

    选取的 k 条测量（1～12）各以独立概率 p_i/1000 缺失（p_i ∈ [0,1000]，
    0‰/1000‰ 即确定事件），未选取测量始终有效。对每种缺失组合：移除该组合中
    的所选测量，重算目标间隙的最紧上/下界；仅当目标有界（双向可达）且整个
    最紧区间落入合格区间 [target_low, target_high] 时计为通过。组合权重
        ∏_{i∈缺失} p_i · ∏_{i∉缺失} (1000−p_i)
    （公分母 1000^k）；通过权重和、未通过权重和分别约为既约分数返回，二者
    之和恰为 1。

    增量更新：底座（未选测量）先做一次全点 Floyd 并按“终端 = 目标点对 +
    各所选测量端点”汇总为终端图；所选测量逐条化作两条边，沿子集递归树加入
    （单条新边 u→v 的更新 d'[i][j] = min(d[i][j], d[i][u]+w+d[v][j])），
    不进入缺失集的分支共享父矩阵、不复制，进入的分支复制后更新，故不为
    2^k 种组合各跑一次全点最短路。
    """
    n = len(points)
    idx = {name: i for i, name in enumerate(points)}
    p = idx[target_from]
    q = idx[target_to]
    sel_ids = {it["id"] for it in items}
    base_constraints = [c for c in constraints if c["id"] not in sel_ids]
    sel_constraints = items  # 已按录入顺序排列

    base: dict = {
        "consistent": True,
        "points": points,
        "targetFrom": target_from,
        "targetTo": target_to,
        "targetLow": target_low,
        "targetHigh": target_high,
    }

    # 底座（未选测量，始终有效）全源最短路
    dist0 = all_pairs_shortest(n, build_edges(points, base_constraints))

    # 终端：目标点对 + 各所选测量端点（去重，保持出现顺序）
    terminals = [p, q]
    seen_t = {p, q}
    cons_of = {c["id"]: c for c in constraints}
    for it in sel_constraints:
        c = cons_of[it["id"]]
        for name in (c["from"], c["to"]):
            t = idx[name]
            if t not in seen_t:
                seen_t.add(t)
                terminals.append(t)
    loc = {t: i for i, t in enumerate(terminals)}
    p_loc, q_loc = loc[p], loc[q]

    base_mat = [[dist0[u][v] for v in terminals] for u in terminals]
    extra: list[list[tuple[int, int, int]]] = []
    for it in sel_constraints:
        c = cons_of[it["id"]]
        u = loc[idx[c["from"]]]
        v = loc[idx[c["to"]]]
        extra.append([(u, v, c["high"]), (v, u, -c["low"])])

    k = len(items)
    miss = [it["missingPermille"] for it in items]
    valid = [1000 - m for m in miss]
    denominator = 1000 ** k

    full_bounds: dict = {}
    pass_weight = 0
    fail_weight = 0
    pass_count = 0
    fail_count = 0

    def dfs(b: int, d: list[list[int | None]], weight: int, mask: int) -> None:
        """处理所选测量 b..k-1；d 为终端图最短路（有效测量＝底座 + b 之前不在 mask 者）。

        mask 的第 i 位为 1 表示测量 i 缺失。weight 为到达该状态已确定因素的
        组合权重乘积；叶节点即一种缺失组合的完整权重（分母 1000^k）。
        “缺失”分支不加入新边、共享父矩阵（不复制）；“有效”分支复制矩阵后
        增量加入该测量的两条边。全程遍历全部 2^k 个叶节点（0‰/1000‰ 的
        零权重分支也计数），结构组合数恒为 2^k。
        """
        nonlocal pass_weight, fail_weight, pass_count, fail_count
        if b == k:
            up = d[p_loc][q_loc]
            lo = None if d[q_loc][p_loc] is None else -d[q_loc][p_loc]
            if mask == 0:
                full_bounds["low"] = lo
                full_bounds["upper"] = up
            ok = (up is not None and lo is not None
                  and target_low <= lo and up <= target_high)
            if ok:
                pass_weight += weight
                pass_count += 1
            else:
                fail_weight += weight
                fail_count += 1
            return
        # 测量 b 缺失：不加入其边，权重乘 p_b（0 时该分支概率为 0，仍计数）
        dfs(b + 1, d, weight * miss[b], mask | (1 << b))
        # 测量 b 有效：复制矩阵后增量加入其两条边，权重乘 (1000−p_b)
        d2 = [row[:] for row in d]
        for u, v, w in extra[b]:
            _relax_terminal_edge(d2, u, v, w)
        dfs(b + 1, d2, weight * valid[b], mask)

    dfs(0, [row[:] for row in base_mat], 1, 0)

    if pass_weight + fail_weight != denominator:
        raise RuntimeError("缺失概率不变式被破坏：通过/未通过权重和不等于 1000^k")

    pass_fr = Fraction(pass_weight, denominator)
    fail_fr = Fraction(fail_weight, denominator)
    full_low = full_bounds["low"]
    full_upper = full_bounds["upper"]

    base.update({
        "selectedItems": [
            {"id": it["id"], "entryIndex": it["entryIndex"],
             "missingPermille": it["missingPermille"],
             "validPermille": 1000 - it["missingPermille"]}
            for it in items
        ],
        "itemCount": k,
        "combinationCount": 1 << k,
        "passCount": pass_count,
        "failCount": fail_count,
        "fullBounded": full_low is not None and full_upper is not None,
        "fullLow": full_low,
        "fullHigh": full_upper,
        "passProbability": _fraction_payload(pass_fr),
        "failProbability": _fraction_payload(fail_fr),
    })
    return base


# ============================ 极值共现工况 ============================

MIN_EXTREME_TARGETS = 1
MAX_EXTREME_TARGETS = 10


def validate_extreme_targets(points: list[str], raw_targets) -> list[dict]:
    """校验极值共现目标清单，返回规范化目标列表；失败抛 ValueError(消息列表)。

    每项：{"from", "to", "side"}；from/to 为两个不同的已录入基准点，
    side ∈ {"lower", "upper"}（最紧下界 / 最紧上界）；同一有序点对、同一
    界别不得重复；共 1～10 项。
    """
    if not isinstance(raw_targets, list):
        raise ValueError(["targets 必须为数组：元素为 {from, to, side} 对象"])
    if not (MIN_EXTREME_TARGETS <= len(raw_targets) <= MAX_EXTREME_TARGETS):
        raise ValueError([f"目标数量须为 {MIN_EXTREME_TARGETS}～{MAX_EXTREME_TARGETS} 个"
                          f"（当前 {len(raw_targets)}）"])

    point_set = set(points)
    errors: list[str] = []
    targets: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for i, raw in enumerate(raw_targets):
        prefix = f"目标 {i + 1}"
        if not isinstance(raw, dict):
            errors.append(f"{prefix}：格式不是对象")
            continue
        p_from = raw.get("from")
        p_to = raw.get("to")
        p_from = p_from.strip() if isinstance(p_from, str) else ""
        p_to = p_to.strip() if isinstance(p_to, str) else ""
        side = raw.get("side")
        side = side.strip() if isinstance(side, str) else ""
        if not p_from:
            errors.append(f"{prefix}：起点为空")
        elif p_from not in point_set:
            errors.append(f"{prefix}：起点 {p_from} 尚未录入")
        if not p_to:
            errors.append(f"{prefix}：终点为空")
        elif p_to not in point_set:
            errors.append(f"{prefix}：终点 {p_to} 尚未录入")
        if p_from and p_to and p_from == p_to:
            errors.append(f"{prefix}：起点与终点必须不同")
        if side not in ("lower", "upper"):
            errors.append(f"{prefix}：界别必须为 lower（最紧下界）或 upper（最紧上界）")
        if p_from and p_to and side in ("lower", "upper"):
            key = (p_from, p_to, side)
            if key in seen:
                errors.append(f"{prefix}：与之前目标重复（同一有序点对、同一界别不得重复）")
            seen.add(key)
        targets.append({"from": p_from, "to": p_to, "side": side})
    if errors:
        raise ValueError(errors)
    return targets


def _component_zeros(points: list[str], constraints: list[dict]) -> list[int]:
    """原始测量图的连通分量零点：返回每个点所属分量中最早录入的基准点下标。"""
    n = len(points)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    idx = {name: i for i, name in enumerate(points)}
    for c in constraints:
        ra, rb = find(idx[c["from"]]), find(idx[c["to"]])
        if ra != rb:
            parent[rb] = ra
    zero = [0] * n
    first: dict[int, int] = {}
    for v in range(n):
        r = find(v)
        if r not in first:
            first[r] = v      # 按下标递增遍历，首次出现即分量内最早录入
        zero[v] = first[r]
    return zero


def extreme_scenario(points: list[str], constraints: list[dict],
                     targets: list[dict]) -> dict:
    """极值共现工况（仅在原清单一致时调用；矛盾检测由调用方先行完成）。

    对完整清单求每个目标有序间隙 x[to] − x[from] 的最紧整数界（上界
    dist[from][to]、下界 −dist[to][from]）；把“达到该界”作为等式条件
    （x[to] − x[from] = 界值，即一对反向等权边）叠加到原测量图上。一个
    目标子集能由同一组坐标同时达到 ⇔ 叠加对应等式边后的差分约束系统仍无
    负权闭环。求目标数最多的可行子集，并列时取输入序号列表字典序最小者：
    按子集大小从大到小、同大小按输入顺序字典序扫描（combinations 天然按
    字典序生成，至多 2^10 个子集），首个可行者即所求。随后给出一组满足
    原测量与入选等式的整数坐标：叠加图的可行势（Bellman-Ford，整数）按
    原始测量图的连通分量整体平移，使每个分量最早录入的基准点为零点
    （分量间无约束，平移不破坏任何不等式）。任一目标两点在原图中不连通
    （该侧界无界）时报明目标序号，不生成工况。
    """
    n = len(points)
    idx = {name: i for i, name in enumerate(points)}
    dist = all_pairs_shortest(n, build_edges(points, constraints))

    # 各目标的最紧整数界（完整清单）
    infos: list[dict] = []
    unbounded: list[int] = []
    for i, t in enumerate(targets):
        p = idx[t["from"]]
        q = idx[t["to"]]
        if t["side"] == "upper":
            bound = dist[p][q]                       # 最紧上界
        else:
            d_qp = dist[q][p]
            bound = None if d_qp is None else -d_qp  # 最紧下界
        infos.append({"index": i + 1, "from": t["from"], "to": t["to"],
                      "side": t["side"], "bound": bound})
        if bound is None:
            unbounded.append(i + 1)

    base: dict = {"consistent": True, "points": points, "targets": infos}
    if unbounded:
        # 目标界无界：报明目标序号，不生成工况
        base["bounded"] = False
        base["unboundedIndices"] = unbounded
        return base

    m = len(targets)
    base_edges = build_edges(points, constraints)
    # 每个目标的等式条件：x[to] − x[from] = 界值（两条反向等权边）
    eq_edges: list[tuple[Edge, Edge]] = []
    for t, info in zip(targets, infos):
        p = idx[t["from"]]
        q = idx[t["to"]]
        b = info["bound"]
        eq_edges.append((Edge(p, q, b, "", True), Edge(q, p, -b, "", False)))

    def feasible(sel: tuple[int, ...]) -> bool:
        edges = list(base_edges)
        for i in sel:
            edges.extend(eq_edges[i])
        return find_negative_cycle(n, edges) is None

    # 目标数最多的可行子集；并列时按输入序号列表字典序最小
    # （combinations 按输入顺序字典序生成，同大小首个可行者即所求；
    #  空集必可行——原清单已一致，chosen 必被赋值）
    chosen: tuple[int, ...] = ()
    for k in range(m, -1, -1):
        hit = None
        for combo in combinations(range(m), k):
            if feasible(combo):
                hit = combo
                break
        if hit is not None:
            chosen = hit
            break

    # 工况坐标：原测量 + 入选等式下的一组可行整数坐标；再按原始连通分量
    # 整体平移，使每个分量最早录入的基准点为零点
    edges = list(base_edges)
    for i in chosen:
        edges.extend(eq_edges[i])
    pot = _feasible_potential(n, edges)
    zeros = _component_zeros(points, constraints)
    coords = [pot[v] - pot[zeros[v]] for v in range(n)]

    chosen_set = set(chosen)
    base.update({
        "bounded": True,
        "selectedCount": len(chosen),
        "selectedIndices": [i + 1 for i in chosen],
        "excludedIndices": [i + 1 for i in range(m) if i not in chosen_set],
        "selected": [infos[i] for i in chosen],
        "excluded": [infos[i] for i in range(m) if i not in chosen_set],
        "coordinates": [
            {"point": points[v], "value": coords[v], "zeroPoint": points[zeros[v]]}
            for v in range(n)
        ],
    })
    return base
