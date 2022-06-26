import numpy as np
import conf

from collections import Counter

def locateOp(j,i):
    """
    这是一个定位函数，目的就是为了将工件工序Oij定位到相应的DNA序列长度上面
    因为后面我们要针对每一个Oij和机器Mi创建一个开始时间和结束时间的二维数组
    该二维数组在我们最优插入法解码中起到了关键作用：
    主要是给每个（机器编号，工件工序）记录开始时间和结束时间，方便之后的空闲插入操作
    :param j: 工件编号
    :param i: 工序编号
    :return: DNA序列中的某一个位置
    """
    Ops=-1
    for k in conf.JOBS_OPERATIONS_NUM.keys():
        if j==k:
            return Ops+i
        else:
            Ops+= conf.JOBS_OPERATIONS_NUM[k]
    return Ops

def test_mutate_os(sample_os,sample_ms):
    """
    该工具函数是专门用于检测基于领域搜索变异的OS序列，变异的OS是否还对应MS序列中机器编号
    若不符合该函数会抛出异常（因为相应的（工件，工序）并不适应原有的机器编号，即不在JOBS_MACHINE_CHOICE中）
    所以借助异常处理（KeyError），直接返回False，若不抛异常，则变异是可行的，即可解码的
    将可行的变异最终返回True
    :param sample_os: 要检测是否可行的变异样本OS
    :param sample_ms: 跟随OS的MS序列
    :return: Bool值
    """
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
        try:
            machineNumber = ms_dict[(j2, j_operationNo)]
            #这里采用异常捕获的方式去检测
        except KeyError:
            return False

    #不抛出异常，说明变异样本可行，返回True，由调用部分继续处理。
    return True

def merge_os_ms(os_pop, ms_pop, num=conf.POP_NUM):
    cache = []# 由于内置列表添加元素的简易性，先建一个总缓存，之后再转换成ndarray返回
    ms_index = 0  # MS序列种群的索引，即要与下面遍历的OS种群序列对应
    for os in os_pop:# 遍历OS种群
        cache_list = []  # 设置一个缓存列表，该缓存列表的目的是统计某一个工件出现次数（即对应工序）
        ms_dict = {}  # ms_dict是一个工件工序：选用机器编号的字典
        index = 0
        for j in conf.JOBS_MACHINE_CHOICE.keys():  ##等于遍历每一个工件工序
            ms_dict[j] = ms_pop[ms_index][index]  ##根据MS序列确定该工件工序的机器编号
            index += 1
        ms_index += 1 # 再次循环时确保与os同步
        cache_backup=cache.copy()
        for j in os:  # 开始遍历os序列中的元素，即工件
            cache_list.append(j)
            j_operationNo = Counter(cache_list)[j]  ##Counter函数计数元素在序列中出现次数
            """
            这里用一个缓存列表cache_list遍历添加的原因就是：
            要统计对应的工件出现了几次（即得到工件的工序j_operationNo）
            """
            try:
                machineNumber = ms_dict[(j, j_operationNo)]  ##用对应的工件工序元组作为键，找到对应的机器编号
                cache.append((j, j_operationNo, machineNumber))##完成OS和MS中一个元素的合并
            except KeyError as k:
                cache=cache_backup
                num-=1
                break
        """
        n 是三维数组中种群的数量，sum(dic_jobsOperations.values())是DNA序列长度，
        因为最小单位是三元组所以-1代表的就是长度3
        """
    return np.array(cache).reshape(num, sum(conf.JOBS_OPERATIONS_NUM.values()), 3)

def get_time(popMerged,plot=False):
    """
    :param popMerged: os与ms融合后的三元数组
    :return: 每个os与ms所需要消耗的时间
    """
    min_time=np.inf
    time_fitness=[]#存放每对序列总消耗时间
    machine_nums=len(conf.MACHINE_JOBS_TIME) #机器总数
    j_op_nums=len(conf.JOBS_MACHINE_CHOICE) #所有工件工序总数
    jobs_lastmachine=[0]*len(conf.JOBS_OPERATIONS_NUM)
    #记录该工件上一道工序使用的机器编号

    for k in popMerged:##因为k是种群中的一对DNA（OS，MS），所以遍历all
        """
        下面的开始时间矩阵Start_Time和结束时间矩阵End_Time可以这样理解
        行：机器编号
        列：根据工件工序顺序排列的总数，即O11、O12......Onm
        两个矩阵中的相同位置的元素，就代表该工件工序在该机器编号上（等于坐标定位）的开始和结束时间。
        对这两个时间的记录不难，但是如何根据这两个矩阵确定空闲时间段才是重中之重
        """
        Start_Time = np.zeros(shape=(machine_nums, j_op_nums))#行数等于机器数量，列数等于总工序
        End_Time = np.zeros(shape=(machine_nums, j_op_nums))#行数等于机器数量，列数等于总工序

        for i in k:##同样对工件工序机器编号三元组进行遍历。
            machineN=i[2] #该工件工序所用的机器编号
            jobN=i[0] #工件序号
            job_op=i[1] #工序序号
            machineTime= conf.MACHINE_JOBS_TIME[machineN][(jobN, job_op)]
            ##根据机器编号和工件工序确定消耗时间machinetime
            locateInDna=locateOp(jobN,job_op) ##定位在总工序DNA中的位置
            end_time_mi_list=End_Time[machineN] ##根据机器编号选择结束时间矩阵中一行
            start_time_mi_list=Start_Time[machineN]##根据机器编号选择结束时间矩阵中一行
            M_j=len([i for i in end_time_mi_list if i>0])
            ##M_j是在machineN中的第几道工序
            """
            开始计算每一道工件工序的开始和结束时间，
            第一种情况是，在工件上是第一个工序，在该机器上也是第一个工序，
            所以直接从0的开始时间插入，结束时间就是它消耗的时间，然后结束(下面的条件不执行了)，
            返回继续遍历计算下一对三元组。
            """
            if job_op==1 and M_j==0:
                Start_Time[machineN][locateInDna]=0
                # 在开始时间矩阵中记录本次插入的开始时间
                End_Time[machineN][locateInDna]=machineTime
                # 在结束时间矩阵中记录本次插入的结束时间
                jobs_lastmachine[jobN-1]=machineN
                # 记录该工件本次使用的机器编号

            #第二种情况，不是第一个工序，但是机器上第一个工序
            ##所以要寻找上一道工序的结束时间，作为目前工序的开始时间
            elif job_op>1 and M_j==0:
                #不知道上一道工序在哪个机器就遍历寻找k中的元素（三元组）
                #符合则把机器编号赋值给last_machine_No
                last_machine_No=jobs_lastmachine[jobN-1]
                job_op_last_consume=End_Time[last_machine_No][locateInDna-1]
                #找到工件上一个工序消耗的时间
                Start_Time[machineN][locateInDna]=job_op_last_consume
                End_Time[machineN][locateInDna]=job_op_last_consume+machineTime
                ##找到后更新开始时间和结束时间矩阵！
                jobs_lastmachine[jobN-1]=machineN
            #下面这种情况是最容易遗漏的
            #就是它是工件上的第一个工序，而不是机器上的第一道工序
            #这时候就要在机器上面寻找空闲位置插入了
            #寻找空闲时间，并插入的操作
            #如果空闲时间不满足，直接插到最后
            elif job_op==1 and M_j!=0:
                end_time_list = []  # end_time为机器Mi对应的那行End_Time里不为0的元素list
                end_time_location = []  # 记下不为0的位置
                for k1 in range(len(end_time_mi_list)):
                    a = end_time_mi_list[k1]
                    if a != 0:
                        end_time_list.append(a)
                        end_time_location.append(k1)
                end_time_list.sort()
                #以上步骤就是找到该机器上存在的结束时间并记录，最后再进行排序（从小到大）
                start_time_list = []  # start_time为end_time对应的位置的开始时间
                for i in end_time_location:
                    start_time_list.append(start_time_mi_list[i])
                start_time_list.sort()
                #同样是记录开始时间，再排序（从小到大）
                number = len(end_time_list) ##number有多少忙碌时间段
                busy_time = np.zeros((number, 2), dtype=int)
                ##维度number是忙碌时间段数量，2就是开始时间和结束时间
                for j in range(len(start_time_list)):
                    busy_time[j][0] = start_time_list[j]
                    busy_time[j][1] = end_time_list[j]

                free_time=np.zeros(shape=(number+1,2))
                #生成一个numpy二维数组来存储待选择的可插入空闲时间数组
                for ti in range(len(busy_time)):
                    free_time[ti][1]=busy_time[ti][0]
                    free_time[ti+1][0]=busy_time[ti][1]
                #从上面的遍历可以得知free_time数组的索引0是空闲开始时间，索引1是空闲结束时间

                for fi in range(len(free_time)):
                    #遍历free_time，在其中寻找可以插入的点
                    free_time_len=free_time[fi][1]-free_time[fi][0]
                    #先计算每一段空闲时间的长度free_time_len
                    if free_time_len>=0: #在没有超出空闲时间的范围里寻找插入点
                        if machineTime<=free_time_len:
                            # 如果插入时间小于空闲时间，即可插入！
                            Start_Time[machineN][locateInDna]=free_time[fi][0]
                            #在开始时间矩阵中记录本次插入的开始时间
                            End_Time[machineN][locateInDna]=free_time[fi][0]+machineTime
                            #在结束时间矩阵中记录本次插入的结束时间
                            jobs_lastmachine[jobN - 1] = machineN
                            #记录该工件本次使用的机器编号
                            break #直接结束遍历
                    else:
                        #若空闲时间段内都不满足插入时间要求，则直接插入机器最后
                        Start_Time[machineN][locateInDna]=busy_time[-1][1]
                        #在开始时间矩阵中记录本次插入的开始时间
                        End_Time[machineN][locateInDna]=busy_time[-1][1]+machineTime
                        #在结束时间矩阵中记录本次插入的结束时间
                        jobs_lastmachine[jobN - 1] = machineN
                        #记录该工件本次使用的机器编号
                        break

            else:
                """最普遍的工件工序插入情况"""
                last_machine_No = jobs_lastmachine[jobN-1]
                #寻找上一道工序使用机器编号
                last_time_consume=End_Time[last_machine_No][locateInDna-1]
                #找到上一道工序使用时间
                end_time_list = []  # end_time为机器Mi对应的那行End_Time里不为0的元素list
                end_time_location = []  # 记下不为0的位置
                for k1 in range(len(end_time_mi_list)):
                    a = end_time_mi_list[k1]
                    if a != 0:
                        end_time_list.append(a)
                        end_time_location.append(k1)
                end_time_list.sort()
                # 以上步骤就是找到该机器上存在的结束时间并记录，最后再进行排序（从小到大）
                start_time_list = []  # start_time为end_time对应的位置的开始时间
                for i in end_time_location:
                    start_time_list.append(start_time_mi_list[i])
                start_time_list.sort()
                # 同样是记录开始时间，再排序（从小到大）
                number = len(end_time_list)  ##number有多少忙碌时间段
                busy_time = np.zeros((number, 2), dtype=int)
                ##维度number是忙碌时间段数量，2就是开始时间和结束时间
                for j in range(len(start_time_list)):
                    busy_time[j][0] = start_time_list[j]
                    busy_time[j][1] = end_time_list[j]

                free_time = np.zeros(shape=(number + 1, 2))
                for ti in range(len(busy_time)):
                    free_time[ti][1]=busy_time[ti][0]
                    free_time[ti+1][0]=busy_time[ti][1]

                for fi in range(len(free_time)):
                    # 遍历free_time，在其中寻找可以插入的点
                    free_time_len=free_time[fi][1]-free_time[fi][0]
                    if last_time_consume<free_time[fi][1]:
                    #将上一道工序结束时间与空闲结束时间点进行比较
                        if last_time_consume<=free_time[fi][0]:
                        #如果上一道工序结束时间不大于某一段空闲结束时间，又不大于其开始时间
                            if machineTime <= free_time_len:
                            #再对本次工序消耗时间与空闲时间段进行比较，如果小与，则可以插入
                                Start_Time[machineN][locateInDna] = free_time[fi][0]
                                #在开始时间矩阵中记录本次插入的开始时间
                                End_Time[machineN][locateInDna] = free_time[fi][0] + machineTime
                                #在结束时间矩阵中记录本次插入的结束时间
                                jobs_lastmachine[jobN - 1] = machineN
                                #记录该工件本次使用的机器编号
                                break
                        else:
                        # 如果上一道工序结束时间不大于某一段空闲结束时间，但大于其开始时间
                            if machineTime <= free_time[fi][1]-last_time_consume:
                            #将本次工序消耗时间与（空闲结束时间-上一道工序结束时间）比较
                            #如果结果不大于，则也可以将本次工序直接插入该空闲时间
                                Start_Time[machineN][locateInDna] = last_time_consume
                                #在开始时间矩阵中记录本次插入的开始时间
                                End_Time[machineN][locateInDna] = last_time_consume + machineTime
                                #在结束时间矩阵中记录本次插入的结束时间
                                jobs_lastmachine[jobN - 1] = machineN
                                #记录该工件本次使用的机器编号
                                break

                    elif free_time_len<0:
                        #这是没有可以插入空闲时间的比较结果，可以直接插入机器最后
                        if last_time_consume<=busy_time[-1][1]:
                        #在插入之前还要考虑上一道工序时间的结束时间，如果不大于机器的最后结束时间
                            Start_Time[machineN][locateInDna] = busy_time[-1][1]
                            #则直接把机器最后结束时间记录在开始时间矩阵中
                            End_Time[machineN][locateInDna] = busy_time[-1][1] + machineTime
                            #结束时间也就是机器最后结束时间+本次工序消耗时间
                            jobs_lastmachine[jobN - 1] = machineN
                            #记录该工件本次使用的机器编号
                            break
                        else:
                        #如果大于机器的最后结束时间
                            Start_Time[machineN][locateInDna] = last_time_consume
                            #则直接把上一道工序结束时间记录在开始时间矩阵中
                            End_Time[machineN][locateInDna] = last_time_consume+machineTime
                            #把上一道工序时间+本次工序消耗时间记录在结束时间矩阵中
                            jobs_lastmachine[jobN - 1] = machineN
                            #记录该工件本次使用的机器编号
                            break

        ultimate_consume=np.max(End_Time)
        #计算最大完工时间，Cmax
        total_consume=np.sum(End_Time-Start_Time)
        #计算总消耗时间，即总负荷Wt
        machine_max_consume=0
        #记录任务结束后所有机器中运行时间，即负荷最大的Wm
        for no in End_Time-Start_Time:
            max_tmp=np.sum(no)
            if machine_max_consume<max_tmp:
                machine_max_consume=max_tmp

        time_fitness.append((ultimate_consume,machine_max_consume,total_consume)) #结束时间矩阵中的最大值，就是最后完工时间

        if plot and ultimate_consume<min_time:
            #该部分只在SchedulePlot.py，即画图模块中使用，当参数plot=True时启动，
            #能够增加画图模块所需的全局变量需求
            conf.START_TIME=Start_Time
            conf.END_TIME=End_Time
            min_time=ultimate_consume
    return np.array(time_fitness)
    #返回含有（Cmax，Wm，Wt）的目标值

def get_pareto_fit(X,ms_len):
    test_data=X.copy() #操作过程中需要对形参进行修改，因此先备份
    pareto_rank = [] #该列表存储在同一帕累托层的多目标值三元组
    p_fit = {} #将多目标三元组值作为键，其适应度作为值
    p_level = 0 #计数对应的帕累托层

    while len(test_data) != 0:
        p_level+=1  #进入一次while循环，表明开始进行一个新的帕累托层计算
        index = []  #记录在该层中已经进入帕累托层的个体索引，之后按照索引删除
        for i in range(len(test_data)):
            if len(test_data) == 1:
                pareto_rank.append(test_data[i])
                index.append(i)
                continue

            t1 = test_data[:i, :] #比较时不能与自身比较，所以拼接一个新的numpy数组
            t2 = test_data[i + 1:, :]
            compare_data = np.concatenate([t1, t2])
            condition1 = np.array(compare_data[:, 0] >= test_data[i][0],dtype=np.int16)
            #得到第一个维度目标值的比较结果，将True和False用整数表示，下同
            condition2 = np.array(compare_data[:, 1] >= test_data[i][1],dtype=np.int16)
            #得到第二个维度目标值的比较结果
            condition3=np.array(compare_data[:,2]>=test_data[i][2],dtype=np.int16)
            #得到第三个维度目标值的比较结果
            min_num = np.min(condition1 + condition2+condition3)
            
            #计算比较结果的总和，比计算最小值
            if min_num >=1:
                #如果总和的最小值不下于1(没有一个值劣于种群其他对应目标值)，
                # 则证明，改多目标值不被任何其他个体支配，其是一个非支配解
                pareto_rank.append(test_data[i])
                # 将该多为目标值记录在pareto层中
                index.append(i)
                # 记录索引，因为已经是非支配解，后续将其删除，因为其不能再参与后续非支配排序的比较

        pareto_array = np.array(pareto_rank)
        #下面就是按照论文公式就算每个个体的适应度值。
        for idx in range(len(pareto_array)):
            x1 = pareto_array[idx][0]
            x2 = pareto_array[idx][1]
            x3=pareto_array[idx][2]
            tuple_x1_x2_x3 = (x1, x2,x3)
            if len(pareto_array) == 1:
                p_fit[tuple_x1_x2_x3] = 1

            else:
                pareto_array_c = np.concatenate([pareto_array[:idx, :], pareto_array[idx + 1:, :]])
                pareto_d_arr = np.abs(pareto_array_c[:, 0] - x1) + \
                               np.abs(pareto_array_c[:, 1] - x2)+\
                               np.abs(pareto_array_c[:, 2]-x3) #论文公式(4)
                C_i=np.min(pareto_d_arr) #论文公式(5)
                fit=np.exp(-(p_level-1+np.exp(-C_i/ms_len))/ms_len) #论文公式(6)

                p_fit[tuple_x1_x2_x3] = fit

        pareto_rank = []
        test_data = np.delete(test_data, index, axis=0)

    return p_fit #最终返回一个多维目标值和其相应适应度的字典

