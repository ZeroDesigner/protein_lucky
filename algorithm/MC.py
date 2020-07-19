# -*-coding:utf-8 -*-
#目标求解sin(x)最大值
import math
import numpy as np
#初始化生成二进制的解
def species_origin_np(population_size,chromosome_length=6):
    population=[list(np.random.choice(['0','1'],chromosome_length, p=[0.5,0.5])) for i in range(population_size)]
    return population
def bin(chromosome_list:'list'):
    tmp=''.join(chromosome_list)
    return int(tmp)
    #从二进制到十进制
def translation(population:'list'):
    tmp=[int(str(i),2) for i in population]
    return tmp
    
def function(population:'list',x_max_value:'int'):

    population=[bin(i) for i in population]
    # 暂存种群中的所有的染色体(十进制)
    temporary=np.array(translation(population))
    #近似归一化处理
    min = np.amin(temporary)
    max = np.amax(temporary)    
    function_arr =  (temporary - min)/x_max_value
    function_sin=np.sin(function_arr)
    #sin(x)作为目标函数
    return function_sin
    
if __name__ == '__main__':
    population=species_origin_np(1000)
    function_sin=function(population,3.14)
    print(function_sin.max())  
    print('done')


