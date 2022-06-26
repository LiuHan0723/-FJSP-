
from conf import *
import numpy as np
from collections import Counter
from tools import locateOp

class PopGenerate:
    """
    本类是专门负责生成初始化种群的类
    """
    def __init__(self,GS=0.6,LS=0.3,RS=0.1):
        self.pop_num=POP_NUM #固定种群数量，在conf中可更改
        self.__os_length=sum(JOBS_OPERATIONS_NUM.values()) #基因序列的总长度
        self.os_pop=[] #用于生成OS种群
        self.ms_pop=[] #用于生成MS种群
        self.gs=GS
        self.ls=LS
        self.rs=RS

    def __generate_random_os(self):
        """
        本方法主要是随机生成一条符合工件出现次数的基因，因为工件出现次数决定了其有多少工序！
        :return:一条OS 工件工序基因
        """
        os_list = []  ##需要维护的OS序列
        """
        随机生成OS序列的难点在于，如何确保每个工件在序列中的出现次数等于它对应的工序数量
        """
        job = 1  # 从第一个工件开始
        for i in JOBS_OPERATIONS_NUM.values():  # i是对应工件的工序数量
            os_list.extend([job] * i)  # OS序列先扩充符合工序数量的工件编号
            job += 1  # 工件编号递增，去扩充其相应工序数量在OS序列中
        """
        上面扩展生成的OS序列是顺序排列，没有随机性，
        下面将根据其长度，打乱索引，根据新的索引排列，
        要用到numpy的random中的permutation函数
        """
        index = np.random.permutation(self.__os_length)  ##生成一列新的索引
        tmp_arr = np.array(os_list)[index]  ##根据新索引重新随机排列（打乱）
        return tuple(tmp_arr)  ##必须转换成一个元组！
        # 这样才可哈希（具有唯一性，不可修改性质），我们才能把它放到集合中。

    def __generate_os_pop(self):
        """
        该函数多次利用self.__generate_random_os()函数批量生成n个OS序列，即最终的OS_POP种群
        :return: OS_POP
        """
        s = set()  # 先定义一个集合对象
        while len(s) < self.pop_num:
            s.add(self.__generate_random_os())
        """
        集合的特点是只能存放不同元素，防止generate_random_OS函数生成相同的OS序列，破坏了种群多样性
        所以n定义是set的长度，当满足set的长度n时才是确定生成了n个不同的OS序列
        """
        return list(s)  ##最后将集合返回成列表，列表比集合好进行操作

    def __generate_ms_pop_GS(self):
        ms_pop = [] #该列表存放每一个MS序列
        ms_list = [0] * self.__os_length  #生成一定长度的MS序列
        machine_time_list = np.zeros(len(MACHINE_JOBS_TIME))  # 时间数组
        job_machine_shuffle = list(JOBS_MACHINE_CHOICE.items())

        for num in range(int(self.gs*self.pop_num)): #生成对应GS比例的MS_POP种群规莫
            np.random.shuffle(job_machine_shuffle) #打乱工件工序
            for i, j in job_machine_shuffle:     #i=工件工序，j=机器耗时
                temp = machine_time_list.copy()
                for k in j:
                    temp[k] += MACHINE_JOBS_TIME[k][i]
                j_sorted = sorted(j)        #将机器编号按顺序排列
                min_index = j_sorted[temp[j_sorted].argmin()]  #找到全局耗时最短的机器编号

                for index in range(len(machine_time_list)):
                    if (index != min_index):  #将不是耗时最短的机器编号的时间数组置为0
                        temp[index] = 0
                ms_list[locateOp(*i)]=min_index #根据工件工序定位基因所在位置，并赋值
                machine_time_list += temp #更新时间数组
            machine_time_list=np.zeros(len(MACHINE_JOBS_TIME)) #完成一个MS后，将时间数组重置为0
            ms_pop.append(ms_list) #在MS_POP中添加MS

        return ms_pop  #最终返回GS生成的MS_POP

    def __generate_ms_pop_LS(self):
        ms_pop = [] #该列表存放每一个MS序列
        ms_list = [] #生成一定长度的MS序列
        machine_time_list = np.zeros(len(MACHINE_JOBS_TIME))  # 时间数组
        count = 0 #计数器，计算是到达第几道工序
        for i, j in JOBS_MACHINE_CHOICE.items():  #i=工件工序，j=机器耗时
            count += 1
            temp = machine_time_list.copy()
            for k in j:
                temp[k] += MACHINE_JOBS_TIME[k][i]
            j_sorted = sorted(j)   #将机器编号按顺序排列
            min_index = j_sorted[temp[j_sorted].argmin()]
            for index in range(len(machine_time_list)):
                if (index != min_index):     #将不是耗时最短的机器编号的时间数组置为0
                    temp[index] = 0
            ms_list.append(min_index) #根据工件工序位置赋值
            machine_time_list += temp #更新时间数组
            if (count == JOBS_OPERATIONS_NUM[i[1]]):
                machine_time_list = np.zeros(len(MACHINE_JOBS_TIME))
                # 完成一个工件的MS后，将时间数组重置为0

        ms_pop.append(ms_list) #在MS_POP中添加MS
        ms_pop=ms_pop*int(self.ls*self.pop_num)
        #加工LS部分生成的MS_POP按比例生成要求规模大小
        return ms_pop

    def __generate_ms_pop_RS(self):
        """
        我们的MS（机器选择序列）是选择根据工序O11、O12...Onm排列
        即第一个为第一个工件第一个工序、第一个工件第二个工序......直到最后一个工件最后一个工序
        这样方便我们对MS之后的交叉和变异操作
        :return:MS_POP
        """
        ms_pop = [] ##需要生成维护的MS_POP序列
        for num in range(int(self.rs*self.pop_num)):
            ms_list=[]
            for i in JOBS_MACHINE_CHOICE.values():  ##根据特定工件工序随机选择可用的机器编号
                ms_list.append(np.random.choice(i))
            ms_pop.append(ms_list)

        return ms_pop

    def __generate_ms_pop(self):
        ms_from_gs=self.__generate_ms_pop_GS()
        ms_from_ls=self.__generate_ms_pop_LS()
        ms_from_rs=self.__generate_ms_pop_RS()

        ms_all=ms_from_ls+ms_from_rs+ms_from_gs

        return ms_all

    def set_os_pop(self):
        """
        生成对象实例后可直接调用的公有方法
        :return: OS_POP
        """
        self.os_pop=self.__generate_os_pop()

    def set_ms_pop(self):
        """
        生成对象实例后可直接调用的公有方法
        :return: MS_POP
        """
        self.ms_pop=self.__generate_ms_pop()

    def get_os_ms(self):
        """
        将OS_POP和MS_POP融合，这种融合肯定不是两个基因序列种群简单的拼接，
        而是生成一种方便我们之后容易计算目标值的组合，每个基因位是一个三元组：（工件，工序，机器）
        :return:一个三维的numpy数组：axis0=种群，axis1=个体，axis2=基因位
        """
        self.set_os_pop()  #先生成OS_POP
        self.set_ms_pop()  #先生成MS_POP
        cache = []  # 由于内置列表添加元素的简易性，先建一个总缓存，之后再转换成ndarray返回
        ms_index = 0  # MS序列种群的索引，即要与下面遍历的OS种群序列对应
        for os in self.os_pop:  # 遍历OS种群
            cache_list = []  # 设置一个缓存列表，该缓存列表的目的是统计某一个工件出现次数（即对应工序）
            ms_dict = {}# ms_dict是一个{工件工序：选用机器编号}的字典
            index = 0
            for j in JOBS_MACHINE_CHOICE.keys():  ##等于遍历每一个工件工序
                ms_dict[j] = self.ms_pop[ms_index][index]  ##根据MS序列确定该工件工序的机器编号
                index += 1
            ms_index += 1  # 再次循环时确保与os同步
            for j in os:  # 开始遍历os序列中的元素，即工件
                cache_list.append(j)
                j_operationNo = Counter(cache_list)[j]  ##Counter函数计数元素在序列中出现次数
                """
                这里用一个缓存列表cache_list遍历添加的原因就是：
                要统计对应的工件出现了几次（即得到工件的工序j_operationNo）
                """
                machineNumber = ms_dict[(j, j_operationNo)]  ##用对应的工件工序元组作为键，找到对应的机器编号
                cache.append((j, j_operationNo, machineNumber))  ##完成OS和MS中一个元素的合并
            """
            n 是三维数组中种群的数量，self.__os__length是DNA序列长度，
            """
        return np.array(cache).reshape(self.pop_num, self.__os_length, -1)

if __name__ == '__main__':
    # 对该类进行测试
    pop=PopGenerate()
    print(pop.get_os_ms())
    MS_POP=pop.ms_pop
    OS_POP=pop.os_pop
    print(MS_POP)
    print(OS_POP)
