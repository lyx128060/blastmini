from collections import defaultdict

def build_kmer_index(db_records, k=8):
    """构建 k-mer 哈希索引"""
    index = defaultdict(list)
    for seq_id, seq in db_records:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            index[kmer].append((seq_id, i))
    return index

def find_seeds(query_seq, kmer_index, k=8):
    """搜索候选种子"""
    seeds = []
    for i in range(len(query_seq) - k + 1):
        kmer = query_seq[i:i+k]
        if kmer in kmer_index:
            for subject_id, subject_pos in kmer_index[kmer]:
                seeds.append((subject_id, subject_pos, i))
    return seeds

def extend_direction(query_seq, subject_seq, q_pos, s_pos, step):
    """向单一方向(左或右)进行无间隙延伸计算"""
    score, max_score = 0, 0
    best_q, best_s = q_pos, s_pos
    matches, curr_matches = 0, 0
    
    curr_q, curr_s = q_pos + step, s_pos + step

    while 0 <= curr_q < len(query_seq) and 0 <= curr_s < len(subject_seq):
        if query_seq[curr_q] == subject_seq[curr_s]:
            score += 2   # 匹配 +2
            curr_matches += 1
        else:
            score -= 1   # 错配 -1

        if score < 0:    # 局部分数低于0或到达边界时停止
            break

        if score > max_score:
            max_score = score
            best_q = curr_q
            best_s = curr_s
            matches = curr_matches

        curr_q += step
        curr_s += step

    return max_score, best_q, best_s, matches

def run_search(query_records, db_records, k=8, top_n=5):
    """串联：建库 -> Seed匹配 -> 延伸 -> 排序"""
    kmer_index = build_kmer_index(db_records, k)
    db_dict = dict(db_records)
    all_results = []

    for q_id, q_seq in query_records:
        seeds = find_seeds(q_seq, kmer_index, k)
        hits = []
        visited = set()

        for s_id, s_pos, q_pos in seeds:
            s_seq = db_dict[s_id]
            
            # 简单去重：同一个粗略区域的seed只延伸一次，提升效率
            region = (s_id, s_pos // 10)
            if region in visited:
                continue
            visited.add(region)

            seed_score = k * 2
            seed_matches = k

            # 向左延伸
            l_score, l_q, l_s, l_matches = extend_direction(q_seq, s_seq, q_pos, s_pos, -1)
            # 向右延伸
            r_score, r_q, r_s, r_matches = extend_direction(q_seq, s_seq, q_pos + k - 1, s_pos + k - 1, 1)

            total_score = seed_score + l_score + r_score
            total_matches = seed_matches + l_matches + r_matches
            align_len = (r_q - l_q + 1)
            identity = total_matches / align_len if align_len > 0 else 0

            hits.append({
                "query_id": q_id,
                "subject_id": s_id,
                "score": total_score,
                "identity": f"{identity:.2%}", # 转换为百分比
                "query_start": l_q + 1,        # 转换为1-based坐标
                "query_end": r_q + 1,
                "subject_start": l_s + 1,
                "subject_end": r_s + 1,
                "pvalue": 0.001                # P-value 占位符，满足 TSV 列格式要求
            })

        # 对当前 query 的结果按最高分排序，截取前 Top N
        hits.sort(key=lambda x: x["score"], reverse=True)
        all_results.extend(hits[:top_n])

    return all_results