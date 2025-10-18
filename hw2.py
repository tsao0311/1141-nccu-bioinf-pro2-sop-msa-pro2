import pandas as pd
import sys

# 增加遞迴深度以防大型測試案例需要大量堆疊
sys.setrecursionlimit(3000)

# -------------------------------------------------------------
# 1. 檔案讀取函數：嚴格 I/O 處理
# -------------------------------------------------------------
def parse_fasta(file_path: str) -> list[str]:
    """ 解析 FASTA 檔案，確保序列是大寫、無空格且長度一致。 """
    sequences = []
    current_sequence = ""
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('>'):
                    # 處理前一條序列
                    if current_sequence:
                        sequences.append(current_sequence.upper().replace(' ', '')) 
                    current_sequence = ""
                else:
                    # 連接序列內容並移除所有空格
                    current_sequence += line.replace(' ', '')
                    
            # 處理最後一條序列
            if current_sequence:
                sequences.append(current_sequence.upper().replace(' ', ''))
                
        if not sequences: return []

        # 檢查序列長度是否一致 (MSA 檔案的必要條件)
        L = len(sequences[0])
        if not all(len(s) == L for s in sequences): 
            # 如果長度不一致，則檔案格式錯誤
            return []
            
    except Exception: 
        # 檔案讀取或處理錯誤
        return []
    
    return sequences

def parse_score_matrix(file_path: str) -> pd.DataFrame:
    """ 解析分數矩陣檔案，返回 DataFrame。 """
    try:
        # 使用 read_csv 處理空格分隔和註釋行
        score_df = pd.read_csv(
            file_path,
            sep=r'\s+',
            comment='#',
            index_col=0,
            skipinitialspace=True,
            engine='python' # 使用 python 引擎處理複雜分隔符
        )
        # 確保所有數據都是整數
        score_df = score_df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    except Exception: 
        return None
    return score_df

# -------------------------------------------------------------
# 2. 核心計算函數：獨立追蹤仿射間隙 + 雙間隙懲罰 0 (V16)
# -------------------------------------------------------------
def calculate_SoP(input_path: str, score_path: str, gopen: int, gextend: int) -> int:
    """
    計算多重序列比對 (MSA) 的 Sum-of-Pair (SoP) 分數。
    採用四種情況分離的獨立追蹤仿射間隙模型，雙間隙 ('-' vs '-') 懲罰為 0。
    """
    sequences = parse_fasta(input_path)
    score_matrix = parse_score_matrix(score_path)

    if not sequences or score_matrix is None:
        return 0

    num_sequences = len(sequences)
    alignment_length = len(sequences[0])
    total_score = 0

    # 遍歷所有不重複的序列對 (i, j)
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            seq1 = sequences[i]
            seq2 = sequences[j]
            part_score = 0

            # 獨立追蹤序列 1 和序列 2 的間隙狀態
            # in_gap1: 追蹤序列 2 插入（序列 1 有間隙）的狀態
            # in_gap2: 追蹤序列 1 插入（序列 2 有間隙）的狀態
            in_gap1 = False  
            in_gap2 = False  

            for k in range(alignment_length):
                char1 = seq1[k]
                char2 = seq2[k]
                score = 0

                # (1) 兩個都是胺基酸 (AA vs AA)
                if char1 != '-' and char2 != '-':
                    try:
                        score += score_matrix.loc[char1, char2]
                    except KeyError:
                        print("Waring")
                    in_gap1 = False 
                    in_gap2 = False
                
                # (2) 序列 1 有間隙，序列 2 是胺基酸 (- vs A)
                elif char1 == '-' and char2 != '-':
                    if in_gap1:
                        score += gextend
                    else:
                        score += gopen
                        
                    in_gap1 = True
                    in_gap2 = False
                
                # (3) 序列 1 是胺基酸，序列 2 有間隙 (A vs -)
                elif char1 != '-' and char2 == '-':
                    if in_gap2:
                        score += gextend
                    else:
                        score += gopen
                        
                    in_gap2 = True
                    in_gap1 = False
                
                # (4) 兩個都是間隙 ('-' vs '-')
                else: 
                    if (in_gap1 or in_gap2):
                        score += gextend
                    else:
                        score += gopen

                    in_gap1 = True
                    in_gap2 = True
                part_score += score
    total_score += part_score
    print(total_score)
    return int(total_score)

if __name__ == "__main__":
    calculate_SoP("examples/test1.fasta", "examples/pam250.txt", -10, -2)
    calculate_SoP("examples/test2.fasta", "examples/pam100.txt", -8, -2)
    # print(calculate_SoP("examples/test1.fasta", "examples/pam250.txt", -10, -2))
    # print(calculate_SoP("examples/test2.fasta", "examples/pam100.txt", -8, -2))
