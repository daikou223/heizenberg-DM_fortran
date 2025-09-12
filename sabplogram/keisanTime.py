import math
time = 0
for i in range(1,15):
    time += math.comb(30,i)/math.comb(27,10)*100/60
print(time)