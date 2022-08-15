num = []
fp = open('./HW4.1.txt', "r", encoding="utf-8-sig")  # 讀取test檔案
# print(lines)
for i in fp.readlines():
    temp = i.strip('\n')
    # print(temp)
    temp = "".join(temp)
    num.append(int(temp))
# print(num)
num.sort(reverse=True)
total_length = sum(num)
N50_sum = 0
index = -1
while(N50_sum <= total_length/2):
    index += 1
    N50_sum += num[index]
print(total_length, total_length/2, num[index])
with open('./HW4.1_answer.txt', 'w') as output:
    output.write(f"Total contig lengh :{total_length}\n")
    output.write(f"50 % total contig lengh::{total_length/2}\n")
    output.write(f"N50 = {num[index]}")
