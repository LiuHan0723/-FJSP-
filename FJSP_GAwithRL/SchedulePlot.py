import matplotlib.pyplot as plt

import numpy as np
import conf
from tools import merge_os_ms,get_time,locateOp
import re

def plot_schedule(os_pop,ms_pop):

    merge_pop=merge_os_ms(os_pop,ms_pop)

    pop_times=get_time(merge_pop,plot=True)[:,0]

    min_idx=np.argmin(pop_times)

    machine_jobs_matrix=np.zeros(shape=(len(conf.MACHINE_JOBS_TIME), len(conf.JOBS_MACHINE_CHOICE)))

    min_pop=merge_pop[min_idx,:,:]
    if(len(conf.JOBS_OPERATIONS_NUM)/len(conf.MACHINE_JOBS_TIME)>=2):
        plt.figure(figsize=(48,12))
    else:
        plt.figure(figsize=(16,8))

    W_total=0
    W_max_machine=[0]*len(conf.MACHINE_JOBS_TIME)

    for m in min_pop:
        m_no=m[-1]
        j,i=m[:2]
        location=locateOp(j,i)
        machine_jobs_matrix[m_no,location]= conf.MACHINE_JOBS_TIME.get(m_no).get((j, i))
        time=int(machine_jobs_matrix[m_no,location])
        W_total+=time
        W_max_machine[m_no]+=time
        plt.barh(m_no, time,
                 left=conf.START_TIME[m_no, location],
                 color="white", edgecolor="black", height=0.5)
        plt.text(conf.START_TIME[m_no, location] + 0.1, m_no, str(j) + "-" + str(i))
        plt.text(conf.START_TIME[m_no, location] + 0.2, m_no - 0.15, str(time),
                 color="red", fontweight="semibold")
        plt.yticks(list(conf.MACHINE_JOBS_TIME.keys()))
        plt.xlabel("Time",size=20,color="red")
        plt.ylabel("MachineNo",size=20,rotation=90)
        i+=1

    machine_num=max(conf.MACHINE_JOBS_TIME)
    C_max=np.max(conf.END_TIME)
    W_max=np.max(W_max_machine)
    W_max_y=np.argmax(W_max_machine)
    plt.vlines(x=C_max,ymin=0,ymax=machine_num,color="g")
    plt.text(C_max,machine_num/2,s="C_max:{}".format(C_max),fontweight="semibold")

    plt.vlines(x=W_max,ymin=W_max_y-0.5,ymax=W_max_y+0.5,color="g")
    plt.text(W_max, W_max_y,s="W_max:{}".format(W_max), fontweight="semibold")

    plt.text(C_max,machine_num,s="W_total:{}".format(W_total), fontweight="semibold")
    plt.title(re.findall("[bK].*\.", conf.PATH)[0] + "FJSP-Gantt chart", fontweight="bold")
    plt.show()

if __name__ == '__main__':
    from PopGenerate import PopGenerate
    pop=PopGenerate()
    pop.get_os_ms()
    conf.OS_POP=pop.os_pop
    conf.MS_POP=pop.ms_pop

    plot_schedule(conf.OS_POP, conf.MS_POP)
    print(conf.START_TIME)

    # plt.barh(1,2,left=1)
    # plt.text(2,1,"(1,1):3")
    # plt.barh(2,1,left=0)
    # plt.show()
