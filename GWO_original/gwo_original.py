import numpy
import matplotlib.pyplot as plt

'''
原始GWO算法
根据在GitHub上找到的matlab的代码改的python代码
GWO论文原文链接：https://www.sciencedirect.com/science/article/pii/S0965997813001853?via%3Dihub
natlab代码链接：https://github.com/search?l=MATLAB&q=GWO&type=Repositories，代码大同小异，


'''


class GWO:
    '''
    # __init__作为类的初始化函数，
    # 在主函数实例化类时若未加入参数，则使用默认函数进行初始化
    # 若需要对初始化狼群则在其函数内根据自己规则进行编写
    '''

    def __init__(self, search_agents=15, search_dim=2, search_iterator=100, lb=-16, ub=16):
        self.search_agents = search_agents
        self.search_dim = search_dim
        self.search_iterator = search_iterator
        self.lb = lb
        self.ub = ub

        self.wolfs = numpy.random.random(size=[self.search_agents, self.search_dim]) * (ub - lb) + lb
        self.alpha_wolf = self.beta_wolf = self.delta_wolf = numpy.zeros(self.search_dim)
        self.alpha_score = self.beta_score = self.delta_score = numpy.inf

    '''
    # 狼群初始化函数，
    # 若要对狼群初始化进行改进可在该函数内编写，
    # 并在__init__函数调用，
    # 或是在hunting函数开头调用
    '''
    def init_wolfs(self, wolf_type=None):
        if wolf_type is None:
            print('initialize wolfs with default type（random）')
        if wolf_type == 'type1':
            print('please code the rule to initialize wolfs in this switch')

    '''
    # 适应度计算函数，
    # 针对每个回合，狼群中每一个个体计算当前适应度，返回适应值进行排序,
    # 主要问题是适应度方程的确定与编写,
    # 该函数中适应度方程为Matlab上默认F10函数
    # 后续方程根据自己理解需要自行编写
    '''
    def fitness(self, each_wolf=None, F=None, prey=None):
        fitness = numpy.inf

        # function1
        if F is None or F == 'F10':
            y = 20 * (1 - numpy.exp(-0.2 * numpy.sqrt(numpy.mean(each_wolf * each_wolf))))
            z = numpy.exp(1) * (1 - numpy.exp(numpy.mean(numpy.cos(2 * numpy.pi * each_wolf)) - 1))
            fitness = y + z

        # function2
        if F == 'F1':
            print(1)

        '''
        # function3

        # function4

        # function5

        '''
        return fitness

    '''
    # 想了一下，按道理选取最优值和位置更新规则一般都是默认由算法定义规则决定，
    # 但是感觉位置更新规则也能进行一定改动，因此将其写入单独函数change_position并在hunting函数中调用
    '''
    def change_position(self, iterator_i):

        a = 2 * (1 - (iterator_i / self.search_iterator))

        for wolf_i in range(self.search_agents):
            r1 = numpy.random.random(size=self.search_dim)
            r2 = numpy.random.random(size=self.search_dim)
            A1 = numpy.array(2 * a * r1 - a)
            C1 = numpy.array(2 * r2)
            alpha_dis = numpy.abs(numpy.multiply(C1, self.alpha_wolf) - self.wolfs[wolf_i])
            X1 = self.alpha_wolf - A1 * alpha_dis

            r1 = numpy.random.random(size=self.search_dim)
            r2 = numpy.random.random(size=self.search_dim)
            A2 = numpy.array(2 * a * r1 - a)
            C2 = numpy.array(2 * r2)
            beta_dis = numpy.abs(numpy.multiply(C2, self.beta_wolf) - self.wolfs[wolf_i])
            X2 = self.beta_wolf - A2 * beta_dis

            r1 = numpy.random.random(size=self.search_dim)
            r2 = numpy.random.random(size=self.search_dim)
            A3 = numpy.array(2 * a * r1 - a)
            C3 = numpy.array(2 * r2)
            delta_dis = numpy.abs(numpy.multiply(C3, self.delta_wolf) - self.wolfs[wolf_i])
            X3 = self.delta_wolf - A3 * delta_dis

            self.wolfs[wolf_i] = (X1 + X2 + X3) / 3
            # end wolf_i

    '''
    # 狩猎函数，主要函数，描述狼群狩猎过程，即寻优的过程，
    # 函数主要有一个大迭代循环和一个灰狼判优循环，以及位置更新函数，画图函数
    # 灰狼群体循环会处理越界灰狼并根据当前狼群适应情况更新领导信息，
    '''
    def hunting(self):

        score_list = []  # 记录每次迭代的最优值

        for iterator_i in range(self.search_iterator):

            for agent_i in range(self.search_agents):

                flag4ub = numpy.array(self.wolfs[agent_i] > self.ub)
                flag4lb = numpy.array(self.wolfs[agent_i] < self.lb)
                self.wolfs[agent_i] = numpy.array(
                    self.wolfs[agent_i] * (~(flag4lb + flag4ub)) + flag4ub * ub + flag4lb * lb)

                fitness = self.fitness(self.wolfs[agent_i])

                if fitness < self.alpha_score:
                    self.alpha_score = fitness
                    self.alpha_wolf = self.wolfs[agent_i]
                elif fitness > self.alpha_score and fitness < self.beta_score:
                    self.beta_score = fitness
                    self.beta_wolf = self.wolfs[agent_i]
                elif fitness > self.alpha_score and fitness > self.beta_score and fitness < self.delta_score:
                    self.delta_score = fitness
                    self.delta_wolf = self.wolfs[agent_i]

            # end agent_i

            self.change_position(iterator_i=iterator_i)

            score_list.append(self.alpha_score)

            print("the ", iterator_i, " times, the best score is: ", self.alpha_score)
            # print('alpha_wolf:',self.alpha_wolf)

        self.plot_score_line(score_list)

        # end iterator_i

    '''
    # 无实际作用，
    # 主要用于在main中输出部分信息，确定GWO类执行过程中部分信息变化情况
    '''
    def show_info(self):
        print(self.fitness(self.wolfs[5]))
        print(self.wolfs)

    '''
    # 画出分数变化曲线，也可自己编写画其他曲线，
    # 在hunting函数末尾调用
    '''
    def plot_score_line(self, score_list):
        x = list(range(len(score_list)))
        plt.plot(x, score_list)
        plt.show()


'''
优化函数之一,别人的
def six_hump_camel_back(variables_values = [0, 0]):
    func_value = 4*variables_values[0]**2 - 2.1*variables_values[0]**4 + (1/3)*variables_values[0]**6 + variables_values[0]*variables_values[1] - 4*variables_values[1]**2 + 4*variables_values[1]**4
    return func_value
'''

if __name__ == '__main__':
    '''
    # 如果已经有数据，读取数据并获取其属性作为搜索维度
    # 狼群数量一般默认10-15，可自行定义
    # 搜索迭代次数自定义
    # 主要关注搜索维度和上下界，需要与数据一致
    # 搜索上下界自定义，很多情况下会将数据归一化从而将上下界定为[-1,1]
    '''
    search_agents = 15
    search_dim = 15
    search_iterator = 500
    lb = -32
    ub = 32

    gwo = GWO(search_agents=search_agents, search_dim=search_dim, search_iterator=search_iterator, lb=lb, ub=ub)
    # gwo.show_info()
    gwo.hunting()

    # x = numpy.random.random(size=5) * 20 - 10




