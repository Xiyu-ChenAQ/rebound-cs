# Cosmic Stars — Physics Engine

> 基于 [REBOUND](https://github.com/hannorein/rebound) 的 Cosmic Stars 专属天体物理计算库。
> 将 REBOUNDx 的物理扩展模块直接内嵌进 REBOUND，去除 Python 依赖，支持 MSVC 编译，面向 C# P/Invoke 分发。

---

## 项目结构

```
rebound-cs/
├── src/            REBOUND 核心 C 源码（N 体积分器）
├── cs/             Cosmic Stars 物理扩展层
│   ├── cs.h                统一对外头文件，用户只需 include 这一个
│   ├── cs_simulation.h/c   框架核心：模块注册、回调调度、生命周期管理
│   ├── cs_gr.h/c           广义相对论后牛顿修正（三种精度模式）
│   └── ...                 （更多模块持续添加中）
├── reboundx/       REBOUNDx 原始 C 源码（仅作参考，不参与编译）
└── examples/       REBOUND C 示例
```

---

## 为什么做这个

REBOUNDx 在 MSVC 上编译失败（VLA、GCC 扩展等问题），导致无法直接集成进
基于 .NET / C# 的 Cosmic Stars 项目。本库的目标是：

- 将 REBOUNDx 的物理模块移植为 **MSVC 兼容的纯 C99 代码**
- 去除所有 Python 绑定和 REBOUNDx 框架层（参数链表、动态注册机制）
- 提供简洁的 C API，方便通过 P/Invoke 从 C# 调用
- 所有物理模块统一通过 `cs/` 层管理，不侵入 REBOUND 核心

---

## 已实现模块

| 模块 | 文件 | 说明 |
|------|------|------|
| GR Potential | `cs/cs_gr.c` | 简单后牛顿势，最快，WHFast 辛安全 |
| GR Single | `cs/cs_gr.c` | 单源完整 1PN 修正（Anderson 1975） |
| GR Full | `cs/cs_gr.c` | 全体天体两两 1PN 修正（Newhall 1983），已修复 VLA |

## 计划中的模块

- `cs_mass` — 恒星质量损失 / 增长
- `cs_radiation` — Poynting-Robertson 辐射拖曳
- `cs_harmonics` — J2 / J4 / J6 引力矩
- `cs_migration` — 轨道迁移（Type I / 力阻尼）
- `cs_tides` — 潮汐（常数时间滞后 / 自旋耦合）

---

## 快速开始

```c
#include "cs/cs.h"

int main(void) {
    // 创建 REBOUND 仿真
    struct reb_simulation* sim = reb_simulation_create();

    // 附加 Cosmic Stars 扩展层
    cs_simulation_t* cs = cs_simulation_create(sim);

    // 开启广义相对论修正（AU + yr + M_sun 单位制）
    cs_enable_gr(cs, CS_GR_POTENTIAL, CS_C_AU_YR_MSUN);

    // 添加天体：恒星 + 行星
    reb_simulation_add_fmt(sim, "m", 1.0);
    reb_simulation_add_fmt(sim, "m a e", 3e-6, 1.0, 0.01);

    // 积分 1000 年
    reb_simulation_integrate(sim, 1000.0);

    // 释放（cs 会随 sim 自动释放）
    reb_simulation_free(sim);
    return 0;
}
```

---

## 上游致谢

本项目建立在两个杰出开源项目的肩膀上。

### REBOUND

**作者：Hanno Rein、Dan Tamayo 及众多贡献者**
https://github.com/hannorein/rebound

REBOUND 是一个极为优雅的 N 体积分器。它的代码结构清晰、接口设计克制、
物理实现严谨，是少数能让人在阅读源码时感到愉悦的科学计算项目之一。
`src/` 目录下的全部核心代码均来自 REBOUND，我们未对其做任何修改。

核心论文：
- Rein & Liu 2012, A&A 537, A128 — 代码结构与引力/碰撞算法
- Rein & Spiegel 2015, MNRAS 446, 1424 — IAS15 高精度积分器
- Rein & Tamayo 2015, MNRAS 452, 376 — WHFast 辛积分器

### REBOUNDx

**作者：Dan Tamayo、Hanno Rein、Pengshuai Shi 及贡献者**
https://github.com/dtamayo/reboundx

REBOUNDx 以极小的侵入性为 REBOUND 提供了丰富的物理扩展模块，
其设计理念（通过回调和 void* 实现零耦合扩展）直接启发了本项目 `cs/` 层的架构。
`reboundx/` 目录保留了其原始 C 源码作为参考。

核心论文：
- Tamayo, Rein, Shi & Hernandez 2020, MNRAS 491, 2885

---

本项目的存在意义不是替代上述项目，而是为了在特定的工程约束下
（MSVC 编译、C# 分发、无 Python 依赖）继续使用它们出色的物理实现。
所有代码遵循 GPL v3 协议。

---

## License

GPL v3 — 详见 [LICENSE](LICENSE)