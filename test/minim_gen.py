from random import choice

n_seq = 1
seq_len = 10

k = 1
win_len = 1 

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

def find_minimizers(seq, k, win_len):
    """ Bruteforce minimizer finding """
    minimizers = []

    for i in range(len(seq) - win_len + 1):
        kmers = []
        for j in range(win_len - k + 1):
            start = i + j
            stop = start + k

            org = seq[start:stop:1]
            comp = get_complement(org)

            kmers.append((convert_kmer(org), start, 0))
            kmers.append((convert_kmer(comp), start, 1))

        minimizers.append(min(kmers))
        kmers.clear()   

    minimizers.sort()

    return minimizers


if __name__ == "__main__":
    # Comment out to override default params
    # n_seq, seq_len, k, win_len = input().split()
    
    seq = gen_seq(seq_len)
    
    print(f"{seq}\n")
    for minim in find_minimizers(seq, k, win_len):
        print(f"{{{minim[0]}, {minim[1]}, {minim[2]}}},")
    
