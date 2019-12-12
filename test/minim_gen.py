from random import choice

seq_len = 15
k = 3
win_len = 5

bases = ('A', 'C', 'G', 'T')

complements = {
    'A': 'T', 'T': 'A',
    'C': 'G', 'G': 'C'
}

masks = {
    'A': '00', 'C': '01',
    'G': '10', 'T': '11'
}

def gen_seq(len):
    return ''.join(choice(bases) for _ in range(len))

def convert_kmer(kmer):
    bin_str = ''.join(masks[base] for base in kmer)
    return int(bin_str, 2)

def get_complement(kmer):
    return ''.join(complements[base] for base in kmer)

def find_minimizers(seq, k, win_len, begin, end):
    """ Bruteforce minimizer finding """

    minimizers = []
    for i in range(begin, end - win_len + 1):
        kmers = []
        for j in range(win_len):
            start = i + j
            stop = start + k

            if stop > len(seq):
                break

            org = seq[start:stop:1]
            comp = get_complement(org)

            kmers.append((convert_kmer(org), start, 0))
            kmers.append((convert_kmer(comp), start, 1))

        if len(kmers) > 0:
            minimizers.append(min(kmers))
        kmers.clear()

    return minimizers

def sorted_and_unique(li):
    ret =  list(set(li))
    ret.sort()

    return ret

if __name__ == "__main__":
    # Comment out to override default params
    seq_len, k, win_len = map(int, input().split())
    
    seq = gen_seq(seq_len)
    minimizers = find_minimizers(seq, k, win_len, 0, seq_len)

    for u in range(1, win_len):
        minimizers.extend(find_minimizers(seq, k, u, 0, u))
    if k < win_len:
        for u in range(k, win_len):
            minimizers.extend(find_minimizers(seq, k, u, seq_len - u, seq_len))

    minimizers = sorted_and_unique(minimizers)
    
    print(f"{seq}\n")
    for minim in minimizers:
        print(f"{{{minim[0]}, {minim[1]}, {minim[2]}}},")
    
