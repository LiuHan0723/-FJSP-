from tools import *
from GeneCrossMutate import GeneCrossMutate


class RL:

    def __init__(self,epislon,gamma,alpha,iternum):
        self.epislon=epislon  #行动策略贪心率
        self.gamma=gamma  #学习折扣率
        self.alpha=alpha  #学习率
        self.iternum=iternum #迭代次数
        self.iternum1=int(iternum/2)  #SARSA算法更新迭代次数
        self.iternum2=int(iternum/2)  #Q-learning算法更新迭代次数
        self.N_STATES = 20  #状态集合数量
        L = np.linspace(0, 1, 21)  #划分状态集[[0,05],[0.05,0.1],[0.1,0.15]...]
        self.STATES = [[L[i], L[i + 1]] for i in range(len(L) - 1)]  #
        PC = np.linspace(0.4, 0.9, 11)  #划分行动集合[[0.4,0.45],[0.45,0.5],[0.5,0.55]...]
        PM = np.linspace(0.01, 0.21, 11) #[[0.01,0.03],[0.03,0.05],[0.05,0.07]...]
        self.ACTIONS = [[PC[i], PC[i + 1], PM[i], PM[i + 1]] for i in range(len(PC) - 1)]
        self.first_fit=[] #收集第一代种群的适应度
        self.q_table=np.zeros(shape=(self.N_STATES,len(self.ACTIONS)))
        #状态-集合数组矩阵的q值表
    def choose_action(self,state):
        q_table=self.q_table
        state_actions = q_table[state,:]
        if (np.random.uniform() > self.epislon) or (state_actions.all() == 0):  # 非贪婪 or 或者这个 state 还没有探索过
            action_no = np.random.choice(range(0,len(self.ACTIONS)))
        else:
            action_no = state_actions.argmax()
            # 贪婪模式
        return action_no

    def set_first_fit(self):

        times=get_time(merge_os_ms(conf.OS_POP, conf.MS_POP)).sum(axis=1)
        # fit_dict = get_pareto_fit(times,len(conf.JOBS_MACHINE_CHOICE))
        max_time = max(times) + 2
        for i in times:
            # self.first_fit.append(fit_dict[tuple(i)])
            self.first_fit.append((max_time - i) / max_time)

    def get_fit(self):

        times = get_time(merge_os_ms(conf.OS_POP, conf.MS_POP)).sum(axis=1)
        # fit_dict = get_pareto_fit(times, len(conf.JOBS_MACHINE_CHOICE))

        max_time = max(times) + 2
        fitness = []  # 需要维护生产的适应度列表

        for i in times:
            # fitness.append(times_fitness_dict[tuple(i)])
            fitness.append((max_time - i) / max_time)

        return fitness

    def get_env_feedback(self,A):
        #行动如何在环境产生交互并最终返回其随之而来的状态和奖励
        last_fit=self.get_fit()
        pc=np.random.uniform(A[0],A[1])
        pm=np.random.uniform(A[2],A[3])
        gene=GeneCrossMutate(p_cross=pc,p_mutate=pm)
        gene.choose()
        gene.mutate()
        conf.OS_POP = gene.next_os_pop
        conf.MS_POP = gene.next_ms_pop
        current_fit=self.get_fit()
        S_=19
        V=sum(current_fit)/sum(self.first_fit)
        for i in range(self.N_STATES):
            if V>=self.STATES[i][0] and V<=self.STATES[i][1]:
                S_=i
                break

        rc=(max(current_fit)-max(last_fit))/max(last_fit)
        rm=(sum(current_fit)-sum(last_fit))/sum(last_fit)
        R=rc+rm

        return S_, R

    def start_learning(self):
        self.set_first_fit()    #得到第一代种群的适应度
        idx_a=np.random.choice(range(len(self.ACTIONS)))  #选择行动集合的索引
        a=self.ACTIONS[idx_a]  #根据索引得到行动集合
        s0, r = self.get_env_feedback(a) #行动与环境交互得到一个新状态s0和奖励r
        idx0=self.ACTIONS.index(a)  #得到该行动在列表中的索引
        for iter1 in range(self.iternum1):
            #SARSA算法更新q值的方式
            idx1=self.choose_action(s0)
            s1,r1=self.get_env_feedback(self.ACTIONS[idx1])


            self.q_table[s0,idx0]=self.q_table[s0,idx0]\
                                 +self.alpha*(r+self.gamma*(self.q_table[s1,idx1])-
                                              self.q_table[s0,idx0])

            s0=s1
            idx0=idx1
            print("第{}次学习完成！".format(iter1+1))
        for iter2 in range(self.iternum2):
            # Q-learning算法更新q值的方式
            idx1 = self.choose_action(s0)
            s1, r1 = self.get_env_feedback(self.ACTIONS[idx1])
            self.q_table[s0, idx0] = self.q_table[s0, idx0] \
                                     + self.alpha * (r + self.gamma * (self.q_table[s1, :].max()) -
                                                     self.q_table[s0, idx0])
            s0 = s1
            idx0 = idx1
            print("第{}次学习完成！".format(self.iternum1+iter2 + 1))



