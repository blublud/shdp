#Three Python packets it needs: numpy,scipy,igraph
#Input file:    a list of PostID that it need to analyze.
#            Like_file: ['PostID','CommentID','UnixSystemTime','LikeFromID','LikeToID']
#Output:    PostID
#        Modularity Value
#        UserID:Degree    for OneGroup
#        UserID:Degree    for AnotherGroup

from sys import stdin, stdout
import numpy as np
from igraph import Graph

from scipy.linalg import eigh  
from scipy.cluster.vq import kmeans2  
from scipy.sparse.linalg import eigen  
from scipy.spatial.kdtree import KDTree  

def rename_clusters(idx):  
    # so that first cluster has index 0  
    num = -1  
    seen = {}  
    newidx = []  
    for id in idx:  
        if id not in seen:  
            num += 1  
            seen[id] = num  
        newidx.append(seen[id])  
    return np.array(newidx)  
    
def cluster_points(L):  
    evals, evcts = eigh(L)  
    evals, evcts = evals.real, evcts.real  
    edict = dict(zip(evals, evcts.transpose()))  
    evals = sorted(edict.keys())  

    for k in range(0,len(evals)):
        if evals[k] > 1e-10:
            startfrom = k
            break
#    print 'startfrom: ', str(evals[k])

    H = np.array([edict[k] for k in evals[startfrom:startfrom+2]]) 
    Y = H.transpose()  
    res, idx = kmeans2(Y, 2, minit='random') 

    return evals[:10000000000], Y, rename_clusters(idx)  


'''
Do the clustering given the like graph.
Input:
    likes: list of like graph's edges
Output:
    ???
'''
def doClustering(likes):

    NodeList = []
    
    for like in likes:
        if like[0] not in NodeList:
            NodeList.append(like[0])
        if like[1] not in NodeList:
            NodeList.append(like[1])

#    print NodeList

    gr = Graph(0)

    gr.add_vertices(len(NodeList))
    #Now, we finish adding nodes.
#    print 'num of node: '
#    print len(NodeList)
    #-------------------------------------
    count = 0
    FromNodeIDtoVerticeID = {}
    for vertice in NodeList:
        FromNodeIDtoVerticeID[vertice] = count
        # In gr.vs, each vertice has two characters.
        gr.vs[count]['verticeID'] = count
        gr.vs[count]['NodeID'] = str(vertice)
        count += 1
#    print 'map from NodeID to verticeID: '
#    print FromNodeIDtoVerticeID 
    #Now we map long NodeID into VerticeID.
    #-------------------------------------
    count = 0
    for line in FromNodeIDtoVerticeID:
        count = count + 1
    #print count
    #-------------------------------------

    for like in likes:
        igraphEdgePair = (FromNodeIDtoVerticeID[like[0]],FromNodeIDtoVerticeID[like[1]])
        gr.add_edges(igraphEdgePair)
    #print summary(gr)

    edgelist = gr.get_edgelist()
    #Now we finish building edges.

    length = len(NodeList)
#    print 'nodelist: ' + str(length)

    #Here we are dealing with special situations in final results.
    if length == 0:
        ClusterResultWithUserID = []
        ClusterResultWithUserID.append('0')

        return ClusterResultWithUserID


    #print '-----------------'
    #Now we are building adjacent matrix for later clustering.
    b = np.arange(0,length*length,1)

    for i in range(0,length*length):
        b[i] = 0
    #print b

    b.shape = length,length

    for i in range(0,len(edgelist)):
        b[edgelist[i][0]][edgelist[i][1]] = b[edgelist[i][0]][edgelist[i][1]] + 1
        b[edgelist[i][1]][edgelist[i][0]] = b[edgelist[i][1]][edgelist[i][0]] + 1
    #Now we finished building adjacent matrix.

    a = [sum(bi) for bi in b]

    G = np.diag(a)  
    L = G - b  
#---------------------------------------------------------------------------------------------------------------------------------
#    for w in range(0,200):
#        evals, Y, idx = cluster_points(L) 
    
#        membership = []

#        for i in range(0,len(idx)):
#            membership.append(str(idx[i]))
#            membership[i] = int(membership[i])

#--------------------------------------------
    checkmodularity = -100
    savemembership = []

    for w in range(0,200):
        evals, Y, idx = cluster_points(L) 
    
        membership = []

        for i in range(0,len(idx)):
            membership.append(str(idx[i]))
            membership[i] = int(membership[i])

        if gr.modularity(membership) > checkmodularity:
            checkmodularity = gr.modularity(membership)
            savemembership = membership


#---------------------------------------------------------------------------------------------------------------------------------
    checkmodularity = gr.modularity(membership)


    OrderIDtoUserID = {}
    #Actually, it is from VerticeIDtoNodeID.
    for line in FromNodeIDtoVerticeID:
        OrderIDtoUserID[FromNodeIDtoVerticeID[line]] = line
    #print OrderIDtoUserID

    #Use OrderIDtoUserID = {}

    OneGroup = {}
    AnotherGroup = {}

    for i in range(0,len(membership)):
        if membership[i] == 0:
            OneGroup[OrderIDtoUserID[i]] = gr.degree(i)

        else:
            AnotherGroup[OrderIDtoUserID[i]] = gr.degree(i)
    
    ClusterResultWithUserID = [OneGroup, AnotherGroup]
    #???Teng also need str(checkmodularity)
    return ClusterResultWithUserID

#end of doClustering
    
#----------------------------------------------------------------------------------------------------------------------------------------

def doClusterToString(likes):
    clusters = ''
    try:            
        ClusterResultWithUserID = doClustering(likes)
        clusters = str(ClusterResultWithUserID)
    except:
        clusters = 'cannot do clustering'
        
    return clusters

        
lastPostId = None
likes = []

for like in stdin:
    [postId, likeFrom, likeTo] = like.strip().split('\t') 
    
    if lastPostId and lastPostId != postId:
        stdout.write(lastPostId)
        clusters = doClusterToString(likes) #pdp: sometimes this takes forever, so print lastPostId first to know which post causes the problem                
        print '\t', clusters        
        likes = []
    
    likes.append((likeFrom, likeTo))
    lastPostId = postId

if lastPostId:
    stdout.write(lastPostId)
    clusters = doClusterToString(likes) #pdp: sometimes this takes forever, so print lastPostId first to know which post causes the problem                
    print '\t', clusters
      
