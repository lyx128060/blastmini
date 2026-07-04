import os

def read_fasta(file_path):
    """
    读取 FASTA 文件并返回生成器 (Generator)。
    使用生成器可以避免将 300MB 的文件一次性全部塞入内存。
    产出: tuple(seq_id, sequence)
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: 找不到文件 {file_path}")

    with open(file_path, 'r') as f:
        seq_id = None
        seq_lines = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq_lines)
                seq_id = line[1:].split()[0] 
                seq_lines = []
            else:
                seq_lines.append(line.upper())
        
        if seq_id is not None:
            yield seq_id, "".join(seq_lines)