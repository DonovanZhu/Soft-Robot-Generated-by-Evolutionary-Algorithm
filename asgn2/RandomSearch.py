import random
import sympy
import numpy

operator_dic = {1: '+', 2: '-', 3: '*', 4: '/', 5: 'sin', 6: 'cos', 7: 'x', 8: 'c'}
operator = ['+', '-', '*', '/']
trangle = ['sin', 'cos']
x_or_const = ['x', 'c']


def generate():
    gene = [0] * 256
    gene[1] = operator_dic[random.randint(1, 4)]
    for i in range(1, 7):
        for j in range(2 ** i, 2 ** (i + 1)):
            if gene[int(j / 2)] in operator:
                gene[j] = operator_dic[random.randint(1, 8)]
            if gene[int(j / 2)] in trangle:
                if j % 2 == 0:
                    gene[j] = operator_dic[random.randint(7, 8)]
    for i in range(2 ** 7, 2 ** 8):
        if gene[int(i / 2)] != 0 and gene[int(i / 2)] not in x_or_const:
            gene[i] = operator_dic[random.randint(7, 8)]
    return gene


def to_string(gene, x):
    result = ""
    if 2 * x > 256 or gene[2 * x] == 0:
        result += gene[x]
    elif gene[x] in trangle:
        result = result + '(' + gene[x] + to_string(gene, 2 * x) + ')'
    else:
        result = result + '(' + to_string(gene, 2 * x) + gene[x] + to_string(gene, 2 * x + 1) + ')'
    return result

