
from parseData import parseBrandimarte,parseKacem
import numpy as np
PATH= "../FJSP_Data/Kacem10x10.txt"
#这里修改数据文件

#针对不同类型数据集进行处理读取
if PATH.startswith("../FJSP_Data/fjsp"):
    JOBS_MACHINE_CHOICE, \
    JOBS_OPERATIONS_NUM, \
    MACHINE_JOBS_TIME = parseBrandimarte(PATH)

if PATH.startswith("../FJSP_Data/Kacem"):
    JOBS_MACHINE_CHOICE, \
    JOBS_OPERATIONS_NUM, \
    MACHINE_JOBS_TIME = parseKacem(PATH)

OS_POP=[]
MS_POP=[]
POP_NUM=100  #这里修改种群规模大小

START_TIME=np.zeros(shape=(len(MACHINE_JOBS_TIME), len(JOBS_MACHINE_CHOICE)))
END_TIME=np.zeros(shape=(len(MACHINE_JOBS_TIME), len(JOBS_MACHINE_CHOICE)))
