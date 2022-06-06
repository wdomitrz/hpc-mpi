import random
import re
import argparse


parser = argparse.ArgumentParser(description='Generate tests.')
parser.add_argument('--genoms', type=int, nargs='+',
                    help='List of lengths of the genoms')
parser.add_argument('--q-num', type=int, help='Number of queries.')
parser.add_argument('--q-len', type=int, help='Length of the queries.')
parser.add_argument('--ensure-query-positive', action='store_true',
                    help='If true, each query has at least one match.')
parser.add_argument('--name', type=str, help='Name of the tests.')
parser.add_argument('--num-tests', type=int, help='Number of tests.')
parser.add_argument('--output', '-o', default=None)


def draw_genom(length, id, test_name='test', letters="ACTG"):
    filename = f'{test_name}_{id}'
    genom = ''.join([random.choice(letters) for _ in range(length)])
    with open(filename, 'w') as f:
        f.write(genom)
    return genom

def draw_query(length, file, query_str=None, letters="ACTG"):
    if query_str is None:
        query = ''.join([random.choice(letters) for _ in range(length)])
    else:
        i = random.randint(0, len(query_str) - length)
        query = query_str[i:i+length]
    file.write(query)
    file.write('\n')
    return query

def add_result(genoms, query, file):
    pattern = re.compile(query)
    result = [str(len(pattern.findall(genom))) for genom in genoms]
    file.write(' '.join(result))
    file.write('\n')
    return result

def maybe_draw(genoms, draw):
    if not draw:
        return None
    return genoms[random.randint(0, len(genoms) - 1)]

if __name__ == '__main__':
    args = parser.parse_args()
    import os

    if args.output is not None:
        dir = f'{args.output}/'
        if not os.path.exists(dir):
            os.makedirs(dir)
    else:
        dir = ''
    for i in range(1, args.num_tests + 1):
        filename = f'{dir}{args.name}{i:03d}'
        genoms = [ draw_genom(size, id, filename)
                    for id, size in enumerate(args.genoms)]
        with open(f'{filename}_queries', 'w') as f:
            queries = [draw_query(args.q_len, f, maybe_draw(genoms,
                                                            args.ensure_query_positive))
                       for _ in range(args.q_num)]
        with open(f'{filename}_results', 'w') as f:
            for q in queries:
                add_result(genoms, q, f)
