质量流计算程序V0.3

包含内容说明：

文本文件：
1. func.py 重要的功能函数
2. main.py 主函数模块
3. mass_flow.txt  main.py 计算得到的结果文件
4. input.txt  输入文件
5. nuclide.txt 所关注核素的卸堆质量

文件夹：
decay
        ---decay.py 48阶CRAM方法堆外衰变计算
				---__init__.py 初始化脚本，用于文件夹外部程序调用
				---example.py 测试例子
				
				--- lib 存放decay.py所需的二进制数据库
				     ---DecayCoeff.npz  各核素的衰变相关系数，包括衰变常数("Lambda")，衰变热系数("DecayHeat")以及放射性毒性系数("Toxicity")
						 ---Matrix_A.npz 衰变系数矩阵

						 
docs 
        --readme.txt  程序文档

figs 存放程序生成的图片

lib 存放主数据库 FCDB.hdf5

test 开发阶段的测试文件夹

其余文件为Python集成开发环境Pycharm 生成的Project文件，可以忽略；

程序目前还需要进一步精简和模块化.
    
