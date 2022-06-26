from PopGenerate import PopGenerate
from tools import *
from SchedulePlot import plot_schedule
from ReinforcementLearning import RL
pop = PopGenerate()  # 生成一个初始化种群的实例对象
pop.get_os_ms()  # 调用get_os_ms()生成OS_POP和MS_POP
# 上面的实例方法最终返回值是OS和MS融合的三元组，生成的OS_POP和MS_POP已经在实例属性中
conf.OS_POP = pop.os_pop  # 将全局配置conf中的OS_POP和MS_POP生成
conf.MS_POP = pop.ms_pop
plot_schedule(conf.OS_POP, conf.MS_POP) #画第一代的甘特图

times_array = get_time(merge_os_ms(conf.OS_POP, conf.MS_POP))
max_idx=times_array.sum(axis=1).argmin()
print("种群初始最优结果(Cmax,Wmax,Wtotal): ",times_array[max_idx])
#第一代种群中适应度最优的子代
print("++++++++++++学习中++++++++++++++")

rl = RL(epislon=0.85, alpha=0.75, gamma=0.2, iternum=50)
#程序主要内容都在这里！更改iternum可减少程序运行时间！
rl.start_learning()  #迭代自学习中

print("==========学习结束==============")
times_array = get_time(merge_os_ms(conf.OS_POP, conf.MS_POP))
#最后一代种群中适应度最优的子代
max_idx=times_array.sum(axis=1).argmin()
print("最后一次学习最优结果(Cmax,Wmax,Wtotal): ",times_array[max_idx])
plot_schedule(conf.OS_POP, conf.MS_POP) #画最后一代的甘特图
