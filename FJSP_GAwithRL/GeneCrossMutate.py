import conf
from tools import merge_os_ms,get_time,test_mutate_os
import numpy as np
from itertools import permutations
from collections import Counter

class GeneCrossMutate:
    """
    该类是在ReinforcementLearning中实例化运行的对象
    """
    def __init__(self,p_cross,p_mutate):
        """
        :param p_cross: 交叉率
        :param p_mutate: 变异率
        """
        self.origin_os_pop= conf.OS_POP #上一代OS种群
        self.origin_ms_pop= conf.MS_POP #上一代MS种群
        self.next_os_pop=[] #需要维护生成的下一代OS种群
        self.next_ms_pop=[] #需要维护生成的下一代MS种群
        self.p_cross=p_cross
        self.p_mutate=p_mutate
        self.length=sum(conf.JOBS_OPERATIONS_NUM.values()) #基因序列长度

    def get_old_merge(self):
        """
        将上一代的种群先融合成方便操作的三元组
        :return: 见tools中的merge_os_ms函数
        """
        return merge_os_ms(self.origin_os_pop,self.origin_ms_pop)

    def get_fitness(self):
        """

        :return: 返回种群的适应度列表
        """
        old_merge=self.get_old_merge() #得到上一代的融合三元组
        times=get_time(old_merge).sum(axis=1) #根据得到的三元组计算每个个体所消耗的时间
        # times_fitness_dict=get_pareto_fit(times,self.length)
        max_time=max(times)+2
        fitness=[] #需要维护生产的适应度列表

        for i in times:
            # fitness.append(times_fitness_dict[tuple(i)])
            fitness.append((max_time-i)/max_time)

        return fitness

    def choose(self):
        """
        该选择过程中包含了交叉过程，中间包括了轮赌盘、和精英保留策略的选择
        """
        fitness=self.get_fitness()
        total_fitness=sum(fitness) #得到所有适应度总和作为基数
        max_fitness_index=np.argmax(fitness) #找到最优秀的基因个体（即适应度最高的）直接进入下一代
        self.next_os_pop.append(list(self.origin_os_pop[max_fitness_index]))
        self.next_ms_pop.append(self.origin_ms_pop[max_fitness_index]) #精英保留策略的体现
        for num in range(conf.POP_NUM - 1):
            """
            由于种群中下一代中已经有了适应度最高的个体，所有还剩POP_NUM-1个基因个体需要生成
            """
            p_c=np.random.random()
              #随机生成一个概率，来决定是否交叉还是直接进入下一代
            #为交叉做准备，根据 个体适应度/种群适应度总和的概率选择两个亲本
            idx = np.random.choice(range(conf.POP_NUM), size=2,
                                   p=np.array(fitness)/total_fitness)
            #若不交叉，则直接继承随机选择的两个亲本中的一个进入下一代
            new_os=np.array(self.origin_os_pop)[idx[0]].tolist()
            new_ms=np.array(self.origin_ms_pop)[idx[0]].tolist()
            #OS满足交叉概率要求：
            if p_c<self.p_cross:
                p1,p2=np.array(self.origin_os_pop)[idx]
                new_os=self.cross_os(p1,p2)
                p1, p2 = np.array(self.origin_ms_pop)[idx]
                new_ms=self.cross_ms(p1,p2)
            #将直接继承的或者交叉之后的个体放入下一代中
            self.next_os_pop.append(new_os)
            self.next_ms_pop.append(new_ms)

    def cross_os(self,p1,p2):
        """
        具体的os交叉步骤
        :param p1: 亲本1
        :param p2: 亲本2
        :return: 交叉后产生的子代OS
        """
        #先随机选择OS序列中的工件，选取数量也随机
        jobs1=np.random.choice(list(conf.JOBS_OPERATIONS_NUM.keys()),
                               size=np.random.randint(1, len(conf.JOBS_OPERATIONS_NUM)),
                               replace=False)

        #以下操作时在两个亲本中定位需要交换工件的索引
        p1_index=np.array([],dtype=np.int16)
        p2_index=np.array([],dtype=np.int16)

        #针对选中的工件遍历其在两个亲代中索引，每次遍历都要拼接
        for job in jobs1:
            jobs1_index=np.argwhere(p1==job).reshape(-1)
            jobs2_index=np.argwhere(p2==job).reshape(-1)
            p1_index=np.concatenate([p1_index,jobs1_index],dtype=np.int16)
            p2_index=np.concatenate([p2_index,jobs2_index],dtype=np.int16)

        p1[list(p1_index)]=p2[list(p2_index)]
        new_os=p1.tolist()

        return new_os

    def cross_ms(self,p1,p2):
        """
        具体的ms交叉步骤
        :param p1: 亲本1
        :param p2: 亲本2
        :return: 交叉后产生的子代MS
        """
        idx=np.random.choice(range(self.length),
                             size=np.random.randint(1,self.length),replace=False)

        p1[list(idx)]=p2[list(idx)]

        new_ms=p1.tolist()

        return new_ms

    def mutate_os(self,sample_os,sample_ms):
        """
        OS采取基于邻域搜索的变异，中间会涉及到全排列，所以产生变异的基因位置不能太多，
        控制在3-6个点。
        :param sample_os: 随机抽取的需要变异的样本OS
        :param sample_ms: 跟随OS的MS基因
        :return:
        """
        idx=np.random.choice(range(self.length),
                             size=np.random.randint(3,6),
                             replace=False
                             )

        origin_array=sample_os[idx]
        mutate_os_pop=[]
        mutate_ms_pop=[]
        ms_list = sample_ms.tolist()
        for i in permutations(origin_array,len(origin_array)):
            sample_os[idx]=i
            sample_os_list=sample_os.tolist()
            if test_mutate_os(sample_os_list,ms_list):
                #这一步就是检测变异的OS基因是否时机器可以加工的工件工序
                #如果是，则通过变异要求，加入mutate_os_pop中
               mutate_os_pop.append(sample_os_list)
               mutate_ms_pop.append(ms_list)

        ##针对在mutate_os_pop中的每一个变异个体，计算他们消耗的时间，选出时间消耗最短的进入下一代
        ##这就是领域搜索的搜索阶段
        mutate_idx = np.argmin(get_time(merge_os_ms(mutate_os_pop, mutate_ms_pop, num=len(mutate_os_pop))))
        new_os = mutate_os_pop[mutate_idx]
        return new_os

    def mutate_ms(self,sample_os,sample_ms):
        """
        这部分代码有一部分和tools中的merge_os_ms相近，但与merge_os_ms有具细微差别
        所以不能直接调用merge_os_ms
        :param sample_os: 随机抽取的需要变异的样本MS
        :param sample_ms: 跟随MS的OS基因
        :return:
        """
        cache=[]
        cache_list = []  # 设置一个缓存列表，该缓存列表的目的是统计某一个工件出现次数（即对应工序）
        ms_dict = {}  # ms_dict是一个工件工序：选用机器编号的字典
        index = 0
        for j1 in conf.JOBS_MACHINE_CHOICE.keys():  ##等于遍历每一个工件工序
            ms_dict[j1] = sample_ms[index]  ##根据MS序列确定该工件工序的机器编号
            index += 1

        for j2 in sample_os:  # 开始遍历os序列中的元素，即工件
            cache_list.append(j2)
            j_operationNo = Counter(cache_list)[j2]  ##Counter函数计数元素在序列中出现次数
            """
            这里用一个缓存列表cache_list遍历添加的原因就是：
            要统计对应的工件出现了几次（即得到工件的工序j_operationNo）
            """
            machineNumber = ms_dict[(j2, j_operationNo)]  ##用对应的工件工序元组作为键，找到对应的机器编号
            cache.append((j2, j_operationNo, machineNumber))  ##完成OS和MS中一个元素的合并
        """
        n 是三维数组中种群的数量，sum(dic_jobsOperations.values())是DNA序列长度，
        因为最小单位是三元组所以-1代表的就是长度3
        """
        one_os_ms=np.array(cache).reshape(1, self.length, -1)

        #上面的操作是融合三元组的部分，先融合成三元组的目的 就是为了能根据（工件，工序）寻找可选的且用时最短的机器
        #下面是具体ms变异操作，先随机生成ms序列的索引
        idx = np.random.choice(range(self.length),
                               size=np.random.randint(1, self.length),
                               replace=False
                               )

        #根据索引，找到需要变异的每个三元组（一个基因位）
        ms_mutate=one_os_ms[0,idx].tolist()
        mutate_machineNo=[] #本列表添加的是所有变异可选机器中时间消耗最短的（贪心算法策略）
        for i in ms_mutate:
            job_op=tuple(i[:2])
            time_list=[]
            for _,d in conf.MACHINE_JOBS_TIME.items():
                    time_list.append(d.get(job_op,999))

            mutate_machineNo.append(np.argmin(time_list))

        #下面的步骤是根据工件工序（1，1）、（1，2）......的顺序还原对应的机器序列顺序
        one_os_ms[0,idx,-1]=mutate_machineNo
        new_ms_list={k:0 for k in conf.JOBS_MACHINE_CHOICE.keys()}
        for k in one_os_ms[0,:,:]:
            key=tuple(k[:2])
            new_ms_list[key]=k[-1]

        new_ms=list(new_ms_list.values())
        return new_ms

    def mutate(self):
        """
        根据变异概率执行的变异操作
        :return: 无返回值。只是将变异后的OS或者MS序列进入下一代self.next_ms_pop中
        """
        for n in range(conf.POP_NUM):
            p_m = np.random.rand()

            if p_m<self.p_mutate:
                sample_os=self.next_os_pop[n]
                sample_ms=self.next_ms_pop[n]
                new_os=self.mutate_os(np.array(sample_os),np.array(sample_ms))
                self.next_os_pop[n]=new_os

                sample_os = self.next_os_pop[n]
                sample_ms = self.next_ms_pop[n]
                new_ms = self.mutate_ms(np.array(sample_os),np.array(sample_ms))
                self.next_ms_pop[n] = new_ms







