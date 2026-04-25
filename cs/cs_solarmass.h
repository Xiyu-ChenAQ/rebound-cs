#ifndef CS_SOLARMASS_H
#define CS_SOLARMASS_H

#include "../src/rebound.h"

#ifdef __cplusplus
extern "C" {
#endif

// 恒星质量损失 / 增长
// 每个时间步后根据粒子的 mdot 参数更新质量
// 挂载在 post_timestep_modifications 回调上
void cs_solarmass(struct reb_simulation* sim);

#ifdef __cplusplus
}
#endif

#endif /* CS_SOLARMASS_H */
