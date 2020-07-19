# -*-coding:utf-8 -*-
#目标求解2*sin(x)+cos(x)最大值
import random
import math
import matplotlib.pyplot as plt
from interval import Interval
import numpy as np
import itertools
#初始化生成chromosome_length大小的population_size个个体的二进制基因型种群
#本质为创造一个列表，其中嵌套列表
#species_origin_np，species_origin_list功能一样
def species_origin_np(population_size,chromosome_length):
    population=[list(np.random.choice(['0','1'],chromosome_length, p=[0.5,0.5])) for i in range(population_size)]
    return population
def species_origin_list(population_size,chromosome_length):
    #第一个列表代表种群列表
    population=[]
    #二维列表，包含染色体和基因
    for i in range(population_size):
        #这个列表代表染色体
        temporary=[]
        #染色体暂存器
        for j in range(chromosome_length):
            temporary.append(str(random.randint(0,1)))
            #随机产生一个染色体,由二进制数组成
        population.append(temporary)
            #将染色体添加到种群中
            #最后返回的是种群
    return population
    
    #将基因拼接起来
    #产生一串二进制数字
def bin(chromosome_list:'list'):
    tmp=''.join(chromosome_list)
    return int(tmp)

    #从二进制到十进制
def translation(population:'list'):
    tmp=[int(str(i),2) for i in population]
    return tmp

    # 目标函数相当于环境 对染色体进行筛选，这里是2*sin(x)+cos(x)
def function(population:'list',max_value:'int'):
    chromosome_length=len(population[0])
    population=[bin(i) for i in population]
    # 暂存种群中的所有的染色体(十进制)
    temporary=translation(population)
    #一个基因代表一个决策变量，其算法是先转化成十进制，然后再除以2的基因个数次方减1(固定值)。
    tmp_x=[temporary[i]*max_value/(math.pow(2,chromosome_length)-1) for i in range(len(temporary))]
    function_list=[math.sin(x) for x in tmp_x]
    #这个x的算法为固定并没有找到解释目前
    #猜想是在做一个类似于标准化的流程
    #这里将2*sin(x)+cos(x)作为目标函数，也是适应度函数
    return function_list    
    
def fitness(function_list): 
    fitness_list=[num if num >= 0 else 0 for num in function_list]
    # 如果适应度小于0,则定为0
    #将适应度添加到列表中
    sum=np.sum(np.array(fitness_list))
    #计算适应度和
    #计算概率
    pro_list=list(map(lambda x: x/sum, fitness_list))
    #计算累积概率
    cum_pro_list=[np.sum(np.array(pro_list[:i+1])) if i > 0 else pro_list[i] for i in range(len(pro_list))]
    #建立区间,下开上闭
    tmp_len=len(cum_pro_list) - 1
    interval_list=[Interval(cum_pro_list[i],cum_pro_list[i+1],upper_closed=True,lower_closed=False) for i in range(tmp_len)]
    interval_list.insert(0,Interval(0,cum_pro_list[0],upper_closed=True,lower_closed=False))
    return fitness_list,interval_list,pro_list
    
    #俄罗斯大轮盘,hahahahahhahh
def roulette(population,interval_list):
    new_population=[]
    #创建一个长度与population的随机列表，元素在(0,1)之间
    dec_pro=list(np.random.random(len(population)))
    #检查随机列表元素所在区间，然后，按照区间选取元素
    for x in range(len(dec_pro)):
        for y in range(len(interval_list)):
            if dec_pro[x] in interval_list[y]:
                new_population.append(population[y])
    return new_population
    #如果一个列表为嵌套列表
    #那么将其展开为一个列表
    #例如：
    # l1=[[1,2],[3,4]]
    # l2=[1,2,3,4]
def to_sim_list(more_list):
    new_list =[] 
    for x in range(len(more_list)):
        for y in range(len(more_list[x])):
            new_list.append(more_list[x][y])
    return new_list
            
    #单点交叉
def crossover(new_population):
    #列表任取两个元素形成新的列表
    tmp_population = list(itertools.combinations(new_population, 2))
    #产生交叉点的随机列表
    tmp_point=[random.randint(0,len(new_population[0])) for i in range(len(tmp_population))]
    #切片交叉
    for i in range(len(tmp_population)):
        tmp_population[i][0][tmp_point[i]:-1],tmp_population[i][1][tmp_point[i]:-1]=tmp_population[i][1][tmp_point[i]:-1],tmp_population[i][0][tmp_point[i]:-1]
    #列表展开
    cross_list = tmp_population = to_sim_list(tmp_population)
    #下面这个可能会报错
    #cross_list =[list(i) for i in list(np.unique(np.array(tmp_population),axis=0))]
    return cross_list
    #变异
    
    
def mutation(cross_list,pm,population_size):
    # 染色体/个体中基因的个数
    for x in range(len(cross_list)):
        for y in range(len(cross_list[0])):
            if random.random()<pm:
                # 如果小于阈值就变异  
                # 生成0到py-1的随机数
                if(cross_list[x][y]=='1'):
                    # 将基因进行单点随机变异，变为0或者1
                    cross_list[x][y]=='0'
                else:
                    cross_list[x][y]=='1'
    #将列表转化为矩阵，然后去重
    cross_list = [list(i) for i in list(np.unique(np.array(cross_list),axis=0))]
    #从列表cross_list中任意选取种群数量大小的元素
    new_population_list=np.random.randint(len(cross_list),size=population_size)
    mutation_list=[cross_list[i] for i in new_population_list]
    new_population=mutation_list
 #  new_population = list(random.sample(mutation_list, 20))
    return new_population
    #迭代
if __name__ == '__main__':
    population=species_origin_np(20,20)
    for i in range(100):
        function_list = function(population,2)
        fitness_list,interval_list,pro_list=fitness(function_list)
        print('This a '+ str(i) +' time:')
        sum=np.sum(np.array(fitness_list))
        print('Score:',sum/20)
        new_population=roulette(population,interval_list)
        cross_list=crossover(new_population)
        new_population=mutation(cross_list,0.01,20)
        population=new_population
    print('done')



