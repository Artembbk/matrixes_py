import copy
from itertools import permutations


class Monom:
    def __init__(self, mult, power):
        self.mult = mult
        self.power = power

    def __str__(self):
        return str(self.mult) + 'x^' + str(self.power)


class MatrixError(BaseException):
    def __init__(self, matrix1, matrix2):
        self.matrix1 = matrix1
        self.matrix2 = matrix2


class Matrix:
    def __init__(self, matrix):
        self.matrix = copy.deepcopy(matrix)

    def __str__(self):
        result = ''
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i]) - 1):
                result += str(self.matrix[i][j]) + '\t'
            result += str(self.matrix[i][-1])
            if i != len(self.matrix) - 1:
                result += '\n'
        return result

    def size(self):
        return len(self.matrix), len(self.matrix[0])

    def __add__(self, other):
        if len(self.matrix) == len(other.matrix) and \
                len(self.matrix[0]) == len(other.matrix[0]):
            result = Matrix([0] * len(self.matrix))
            for i in range(len(self.matrix)):
                result.matrix[i] = [0] * len(self.matrix[0])
            for i in range(len(self.matrix)):
                for j in range(len(self.matrix[i])):
                    result.matrix[i][j] = self.matrix[i][j] + \
                                          other.matrix[i][j]
            return result
        else:
            raise MatrixError(self, other)

    def __mul__(self, other):
        if type(other) == int or type(other) == float:

            result = Matrix([0] * len(self.matrix))
            for i in range(len(self.matrix)):
                result.matrix[i] = [0] * len(self.matrix[0])

            for i in range(len(self.matrix)):
                for j in range(len(self.matrix[i])):
                    result.matrix[i][j] = self.matrix[i][j] * other

            return result

        elif isinstance(other, Matrix) and \
                len(self.matrix[0]) == len(other.matrix):

            result = Matrix([0] * len(self.matrix))
            for i in range(len(self.matrix)):
                result.matrix[i] = [0] * len(other.matrix[0])

            for i in range(len(self.matrix)):
                for j in range(len(other.matrix[0])):
                    for k in range(len(self.matrix[0])):
                        result.matrix[i][j] += self.matrix[i][k] \
                                               * other.matrix[k][j]
            return result
        else:
            raise MatrixError(self, other)

    __rmul__ = __mul__

    def transpose(self):
        tMatrix = Matrix([0] * len(self.matrix[0]))
        for i in range(len(self.matrix[0])):
            tMatrix.matrix[i] = [0] * len(self.matrix)
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                tMatrix.matrix[j][i] = self.matrix[i][j]
        self.matrix = tMatrix.matrix
        return self

    def transposed(self):
        tMatrix = Matrix([0] * len(self.matrix[0]))
        for i in range(len(self.matrix[0])):
            tMatrix.matrix[i] = [0] * len(self.matrix)
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                tMatrix.matrix[j][i] = self.matrix[i][j]
        return tMatrix

    def toLatex(self):
        result = '$\\rightarrow$\n$\\begin{pmatrix}\n'
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i]) - 1):
                result += str(self.matrix[i][j]) + ' & '
            result += str(self.matrix[i][-1])
            if i != len(self.matrix) - 1:
                result += '\\\ \n'
        result += '\n\end{pmatrix}$'
        return result

    def trace(self):
        result = 0
        for i in range(len(self.matrix)):
            result += self.matrix[i][i]
        return result

    def toStep(self, curRaw=0, j=0):
        while self.findNotZeroInColumn(j) == -1:
            j += 1
        if self.matrix[curRaw][j] == 0:
            self.elem2(self.findNotZeroInColumn(j, curRaw + 1), curRaw)
        for i in range(curRaw + 1, len(self.matrix)):
            self.elem1(i, curRaw, -self.matrix[i][j] / self.matrix[curRaw][j])
        if curRaw + 1 < len(self.matrix) and j + 1 < len(self.matrix[0]):
            self.toStep(curRaw + 1, j + 1)

    def elem1(self, i, j, a):
        for k in range(len(self.matrix[0])):
            self.matrix[i][k] += self.matrix[j][k] * a
        return self

    def elem2(self, i, j):
        self.matrix[i], self.matrix[j] = self.matrix[j], self.matrix[i]
        return self

    def elem3(self, i, a):
        for k in range(len(self.matrix[0])):
            self.matrix[i][k] *= a
        return self

    def findNotZeroInColumn(self, j, start=0):
        for i in range(start, len(self.matrix)):
            if self.matrix[i][j] != 0:
                return i
        return -1

    def det(self):
        res = 0
        for permutation in permutations(range(0, self.size()[0])):
            term = sgn(permutation)
            for i in range(len(permutation)):
                term *= self.matrix[i][permutation[i]]
            res += term
        return res

    def detx(self):
        res = []
        for permutation in permutations(range(0, self.size()[0])):
            mult = sgn(permutation)
            power = 0
            for i in range(len(permutation)):
                if type(self.matrix[i][permutation[i]]) == int:
                    mult *= self.matrix[i][permutation[i]]
                else:
                    power += 1
            res.append(Monom(mult, power))
        for i in range(len(res)):
            for j in range(i + 1, len(res)):
                if res[i].power == res[j].power:
                    res[i].mult += res[j].mult
                    res[j].mult = 0
        return res


def readMatrix(f):
    matrix = []
    line = f.readline()
    while line.strip() != '--':
        matrix.append(list(map(int, line.strip().split())))
        line = f.readline()
    return matrix


def toLatexFile(m, f):
    f.write(m.toLatex() + '\n')


def toFile(m, f):
    result = ''
    for i in range(len(m.matrix)):
        for j in range(len(m.matrix[i]) - 1):
            result += str(m.matrix[i][j]) + ' '
        result += str(m.matrix[i][-1])
        if i != len(m.matrix) - 1:
            result += '\n'
    result += '\n----------\n'
    f.write(result)


def sgn(per):
    n = len(per)
    res = 0
    for i in range(n):
        for j in range(i):
            if per[j] > per[i]:
                res += 1
    return (-1) ** (res)


def choose_command():
    global m
    print('Выбери действие:')
    print('-1. Отменить действие (кроме вывода в файлы)')
    print('1. Поменять строчки местами')
    print('2. Умножить строку на скаляр и прибавить к другой строке')
    print('3. Умножить строку на скаляр')
    print('4. Привести к целым числам')
    print('5. Вывести в файлы')
    print('6. Умножить столбец на скаляр')
    print('7. Умножить столбец на скаляр и прибавить к другому столбцу')
    print('8. Поменять столбцы местами')
    print('9. Транспонировать')
    print('10. Выход')
    request = int(input())
    print('ВАЖНО! НУМЕРАЦИЯ С НУЛЯ!')

    if request == -1:
        try:
            previousMatrixes.pop()
            m = copy.deepcopy(previousMatrixes[-1])
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 1:
        try:
            print('Введи номера строк:')
            i, j = map(int, input().split())
            m.elem2(i, j)
            previousMatrixes.append(copy.deepcopy(m))
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 2:
        try:
            print('Строка, к которой надо прибавить')
            i = int(input())
            print('Строка, которую надо прибавить')
            j = int(input())
            print('Скаляр')
            a = int(input())
            m.elem1(i, j, a)
            previousMatrixes.append(copy.deepcopy(m))
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 3:
        try:
            print('Введи номер строки')
            i = int(input())
            print('1. Умножить')
            print('2. Поделить')
            ans = int(input())
            if ans == 1:
                print('На что умножить')
                a = int(input())
                m.elem3(i, a)
            if ans == 2:
                print('На что поделить')
                a = int(input())
                m.elem3(i, 1 / a)
            previousMatrixes.append(copy.deepcopy(m))
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 4:
        for i in range(len(m.matrix)):
            for j in range(len(m.matrix[i])):
                m.matrix[i][j] = int(m.matrix[i][j])
        previousMatrixes.append(copy.deepcopy(m))
    elif request == 5:
        outfileLatex = open("LtxMatrixOut.txt", 'a')
        outfile = open('matrixOut.txt', 'a')
        toLatexFile(m, outfileLatex)
        toFile(m, outfile)
        outfileLatex.close()
        outfile.close()
    elif request == 6:
        try:
            print('Введи номер столбца')
            i = int(input())
            print('1. Умножить')
            print('2. Поделить')
            ans = int(input())
            if ans == 1:
                print('На что умножить')
                a = int(input())
                m.transpose()
                m.elem3(i, a)
                m.transpose()
            if ans == 2:
                print('На что поделить')
                a = int(input())
                m.transpose()
                m.elem3(i, 1 / a)
                m.transpose()
            previousMatrixes.append(copy.deepcopy(m))
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 7:
        try:
            print('Столбец, к которому надо прибавить')
            i = int(input())
            print('Столбец, который надо прибавить')
            j = int(input())
            print('Скаляр')
            a = int(input())
            m.transpose()
            m.elem1(i, j, a)
            m.transpose()
            previousMatrixes.append(copy.deepcopy(m))
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 8:
        try:
            print('Введи номера столбцов:')
            i, j = map(int, input().split())
            m.transpose()
            m.elem2(i, j)
            m.transpose()
            previousMatrixes.append(copy.deepcopy(m))
        except Exception:
            print('Кажется ты что то не то ввел. Попробуй еще раз')
            choose_command()
    elif request == 9:
        m.transpose()
        previousMatrixes.append(copy.deepcopy(m))
    elif request == 10:
        return
    print(m)
    choose_command()


infile = open('matrixIn.txt')
m = Matrix(readMatrix(infile))
print(m)
previousMatrixes = [copy.deepcopy(m)]
# print("Определитель равен:", m.det())
choose_command()
