def parseBrandimarte(path):
    JOBS_MACHINE_CHOICE = {}
    JOBS_OPERATIONS_NUM = {}
    MACHINE_JOBS_TIME = {}
    with open(path, 'r') as f:
        lines = f.readlines()
        head = lines[0].strip().split()

        for i in range(int(head[1])):
            MACHINE_JOBS_TIME[i] = {}  ##生成机器在每个工件工序的时间字典
        job = 0
        for j in lines[1:]:  # 该循环是从第二行开始遍历每一行
            job += 1  # 其中每一行代表一个一个工件
            j_list = j.split("    ")
            op_num = int(j_list[0])
            op = 0
            JOBS_OPERATIONS_NUM[job] = op_num
            for k in j_list[1:]:  # 在每一行（工件）中遍历每一个工序
                op += 1
                k_list = k.split("  ")
                JOBS_MACHINE_CHOICE[(job, op)] = []  # 针对每一个工件工序提前设置字典的KV
                for l in k_list[1:]:  # 遍历每一个工件工序可用的机器编号以及所需时间
                    l_list = l.split()
                    JOBS_MACHINE_CHOICE[(job, op)].append(int(l_list[0]))
                    MACHINE_JOBS_TIME[int(l_list[0])][(job, op)] = int(l_list[1])

    return JOBS_MACHINE_CHOICE, JOBS_OPERATIONS_NUM, MACHINE_JOBS_TIME
    # 最终返回三个设定的全局变量


def parseKacem(path):
    JOBS_MACHINE_CHOICE = {}
    JOBS_OPERATIONS_NUM = {}
    with open(path, "r", encoding="utf8") as f:
        lines = f.readlines()
        machine_num = int(lines[0].split()[1])
        MACHINE_JOBS_TIME = {i: {} for i in range(machine_num)}
        job_no = 0
        for line in lines[1:]:  #该循环是从第二行开始遍历每一行
            job_no += 1        #其中每一行代表一个一个工件
            lin = line.split("  ")
            op_num = int(lin[0])
            op_no = 1         #工序编号
            idx = 0           #读取每一行每个数字的指针
            JOBS_OPERATIONS_NUM[job_no] = op_num
            for i in range(op_num):  #在每一行（工件）中遍历每一个工序
                temp_list = []
                l2 = lin[1].split()
                for j in range(int(l2[idx])):  #遍历每一个工件工序可用的机器编号以及所需时间
                    idx += 1
                    machine_No = int(l2[idx])
                    idx += 1
                    machine_time = int(l2[idx])
                    temp_list.append(machine_No - 1)
                    MACHINE_JOBS_TIME[machine_No - 1][(job_no, op_no)] = machine_time
                JOBS_MACHINE_CHOICE[(job_no, op_no)] = temp_list
                idx += 1 #指针指向下一个工序
                op_no += 1 #工序郑增加
    return JOBS_MACHINE_CHOICE, JOBS_OPERATIONS_NUM, MACHINE_JOBS_TIME
    #最终返回三个设定的全局变量


if __name__ == '__main__':
    # 测试数据
    a, b, c = parseBrandimarte("../FJSP_Data/fjsp.brandimarte.Mk02.m6j10c6.txt")
    print("JOB_MACHINE_CHOICE:")
    print(a)
    print("JOBS_OPERATIONS_NUM:")
    print(b)
    print("MACHINE_JOBS_TIME:")
    print(c)
