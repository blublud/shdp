from sys import stdin

for line in stdin:
    fields = line.strip().split(',')
    postId = fields[0]
    likeFrom = fields[3]
    likeTo = fields[4]
    
    print postId + '\t' + fields[3] + '\t' + fields[4]
    


