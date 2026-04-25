/**
 * @file    cs_radiation.c
 * @brief   辐射压与 Poynting-Robertson 拖曳力 / Radiation pressure and PR drag
 *
 * 基于 Burns, Lamy & Soter (1979) 方程 (5)。
 * 同时处理两个物理效果：
 *   1. 辐射压  — 沿径向向外推，削弱引力，等效为引力乘以 (1 - beta)
 *   2. PR 拖曳 — 切向减速，使粒子轨道缓慢向内螺旋（Poynting-Robertson 效应）
 *
 * Based on Burns, Lamy & Soter (1979), Icarus 40, 1.
 * Ported from REBOUNDx radiation_forces.c by Dan Tamayo & Hanno Rein.
 * Adapted for Cosmic Stars by Xiyu Chen (removed reboundx framework dependency).
 *
 * 速度相关力 — 使用前必须设置：
 *   sim->force_is_velocity_dependent = 1;
 * cs_enable_radiation() 会自动设置此标志。
 */

#include "cs_radiation.h"
#include "cs_simulation.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

/* ============================================================
 * SECTION 1: 底层计算函数
 *
 * 物理推导（Burns et al. 1979 eq. 5）：
 *
 * 辐射源（恒星）发出的光子对尘埃粒子施加两种力：
 *
 * (A) 辐射压（Radiation Pressure）
 *     光子打到粒子上产生的动量转移，方向沿径向向外。
 *     大小 = beta * F_gravity，其中 beta = L*Q_pr/(4*pi*G*M*c*rho*s)
 *     beta 越大，粒子越容易被吹走（beta > 1 时粒子被吹离系统）
 *
 * (B) Poynting-Robertson 拖曳（PR Drag）
 *     运动的粒子看到的光子有多普勒蓝移（前方）和红移（后方），
 *     导致净动量转移有一个切向分量，方向与速度相反，使粒子减速。
 *     这是导致太阳系尘埃盘向内螺旋的主要机制。
 *
 * 合力加速度（在辐射源参考系中）：
 *   a = beta * G * M / r^2 * [(1 - r_dot/c) * r_hat - v_rel/c]
 *
 * 其中：
 *   r_hat  = (r_particle - r_source) / |r_particle - r_source|  径向单位向量
 *   r_dot  = (r_rel · v_rel) / r                                 径向速度（标量）
 *   v_rel  = v_particle - v_source                               相对速度向量
 *   第一项 (1 - r_dot/c) * r_hat  — 辐射压（含一阶多普勒修正）
 *   第二项 -v_rel/c               — PR 拖曳
 * ============================================================ */

void cs_calculate_radiation_forces(
    struct reb_simulation* const sim,
    struct reb_particle* const particles,
    int N,
    double c,
    int source_index)
{
    const struct reb_particle source = particles[source_index];

    /* mu = G * M_source，辐射压大小的基准引力参数 */
    const double mu = sim->G * source.m;

    for (int i = 0; i < N; i++) {
        /* 辐射源自身不受辐射力 */
        if (i == source_index) continue;

        /* 只有设置了 beta 的粒子才受辐射力
         * beta = 0 表示未设置，跳过 */
        cs_particle_params_t* params = cs_particle_params_get(&particles[i]);
        if (params == NULL || params->beta == 0.0) continue;

        const double beta = params->beta;
        const struct reb_particle p = particles[i];

        /* --- 计算相对位置向量（从辐射源指向粒子）--- */
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double dr2 = dx*dx + dy*dy + dz*dz;
        const double dr = sqrt(dr2); /* 粒子到辐射源的距离 */
        if (dr < DBL_EPSILON) continue;

        /* --- 计算相对速度向量 --- */
        const double dvx = p.vx - source.vx;
        const double dvy = p.vy - source.vy;
        const double dvz = p.vz - source.vz;

        /* 径向速度标量 r_dot = (r_rel · v_rel) / r
         * 正值表示粒子正在远离辐射源 */
        const double rdot = (dx*dvx + dy*dvy + dz*dvz) / dr;

        /* 辐射压加速度大小 a_rad = beta * G * M / r^2
         * 这是纯辐射压（无 PR 拖曳修正时）的大小 */
        const double a_rad = beta * mu / (dr * dr);

        /* --- Burns et al. (1979) 方程 (5) ---
         * a_vec = a_rad * [(1 - r_dot/c) * r_hat - v_rel/c]
         *
         * 分解：
         *   (1 - r_dot/c) * r_hat  → 辐射压方向（径向），含多普勒修正
         *   - v_rel / c            → PR 拖曳方向（与相对速度反向）
         *
         * r_hat 的各分量 = (dx/dr, dy/dr, dz/dr) */
        particles[i].ax += a_rad * ((1.0 - rdot/c) * dx/dr - dvx/c);
        particles[i].ay += a_rad * ((1.0 - rdot/c) * dy/dr - dvy/c);
        particles[i].az += a_rad * ((1.0 - rdot/c) * dz/dr - dvz/c);
    }
}

/* ============================================================
 * SECTION 2: 辅助工具 — 计算 beta 值
 *
 * beta 的物理含义：辐射压力与引力之比
 *   beta = F_rad / F_grav
 *        = 3 * L * Q_pr / (16 * pi * G * M * c * rho * s)
 *
 * 其中：
 *   L     — 辐射源光度
 *   Q_pr  — 辐射压效率（完全吸收/散射取 1.0，实际值取决于粒子材质）
 *   G     — 引力常数
 *   M     — 辐射源质量
 *   c     — 光速
 *   rho   — 粒子密度
 *   s     — 粒子半径
 *
 * beta > 1：辐射压 > 引力，粒子被吹离系统
 * beta ~ 0.5：临界逃逸轨道
 * beta << 1：粒子以螺旋方式缓慢向内漂移（PR 拖曳主导）
 * ============================================================ */

double cs_radiation_calc_beta(
    double G,
    double c,
    double source_mass,
    double source_luminosity,
    double radius,
    double density,
    double Q_pr)
{
    return 3.0 * source_luminosity * Q_pr
           / (16.0 * M_PI * G * source_mass * c * density * radius);
}

/* ============================================================
 * SECTION 3: 调度入口 — 由 cs_dispatch_additional_forces() 调用
 *
 * 逻辑：
 *   1. 遍历所有粒子，找到标记了 rad_source = 1 的粒子作为辐射源
 *   2. 如果没有任何粒子标记辐射源，默认使用 particles[0]（通常是恒星）
 *   3. 支持多辐射源（如双星系统），每个辐射源单独调用一次计算
 * ============================================================ */

void cs_radiation_additional_forces(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;

    const int N = sim->N - sim->N_var;
    const double c = cs->rad_c;
    struct reb_particle* const particles = sim->particles;
    int source_found = 0;
    int i;

    /* 寻找辐射源粒子（rad_source == 1）*/
    for (i = 0; i < N; i++) {
        cs_particle_params_t* params = cs_particle_params_get(&particles[i]);
        if (params != NULL && params->rad_source == 1) {
            source_found = 1;
            cs_calculate_radiation_forces(sim, particles, N, c, i);
        }
    }

    /* 没有显式标记辐射源时，默认 particles[0] 为辐射源 */
    if (!source_found) {
        cs_calculate_radiation_forces(sim, particles, N, c, 0);
    }
}