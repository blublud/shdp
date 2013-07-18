#Three Python packets it needs: numpy,scipy,igraph
#Input file:	a list of PostID that it need to analyze.
#	    	Like_file: ['PostID','CommentID','UnixSystemTime','LikeFromID','LikeToID']
#Output:	PostID
#		Modularity Value
#		UserID:Degree	for OneGroup
#		UserID:Degree	for AnotherGroup

import csv
import numpy as np
from igraph import *


import scipy.stats as st
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

#----------------------------------------------------------------------------------------------------------------------------------------

PostList = open('PostListForOccupyLA_PostMoreThan35Comments.txt','r')
rawlikefile = open('occupyla_int_withUnixTime_sorted_addpostID.csv','r');

#rawwritefile = open('To_Sincere_PeopleClusterResultsForEachPost_PostMoreThan35Comments.txt','w');

#PostID = '282370021775789_349187595101578'#about Paul#b
#PostID = '282370021775789_291368864209238'#attitude about LA police#c
#PostID = '282370021775789_187307108026795'#one part#d
#PostID = '282370021775789_384417974904326'#almost two parts

rawlikedata = csv.DictReader(rawlikefile,['PostID','CommentID','UnixSystemTime','LikeFromID','LikeToID'])


#listPostID = ['282370021775789_187307108026795' ,'282370021775789_295392037140254']
listPostID = []
for line in PostList:#in total 1682 posts#400now
	listPostID.append(line.strip())
#print len(listPostID)

ClusterResultWithUserID = {}

for PostID in listPostID:
	PostID = '282370021775789_291368864209238'

	NodeList = []
	rawlikefile.seek(0)
	for rawlikedata_line in rawlikedata:
		if rawlikedata_line['PostID'] == PostID:
			if rawlikedata_line['LikeFromID'] not in NodeList:
				NodeList.append(rawlikedata_line['LikeFromID'])
			if rawlikedata_line['LikeToID'] not in NodeList:
				NodeList.append(rawlikedata_line['LikeToID'])

#	print NodeList

	gr = Graph(0)

	gr.add_vertices(len(NodeList))
	#Now, we finish adding nodes.
#	print 'num of node: '
#	print len(NodeList)
	#-------------------------------------
	count = 0
	FromNodeIDtoVerticeID = {}
	for vertice in NodeList:
		FromNodeIDtoVerticeID[vertice] = count
		# In gr.vs, each vertice has two characters.
		gr.vs[count]['verticeID'] = count
		gr.vs[count]['NodeID'] = str(vertice)
		count += 1
#	print 'map from NodeID to verticeID: '
#	print FromNodeIDtoVerticeID 
	#Now we map long NodeID into VerticeID.
	#-------------------------------------
	count = 0
	for line in FromNodeIDtoVerticeID:
		count = count + 1
	#print count
	#-------------------------------------

	rawlikefile.seek(0)

	for rawlikedata_line in rawlikedata:
		if rawlikedata_line['PostID'] == PostID:
			igraphEdgePair = (FromNodeIDtoVerticeID[rawlikedata_line['LikeFromID']],FromNodeIDtoVerticeID[rawlikedata_line['LikeToID']])
			gr.add_edges(igraphEdgePair)
	#print summary(gr)

	edgelist = gr.get_edgelist()
	#Now we finish building edges.

	length = len(NodeList)
#	print 'nodelist: ' + str(length)

	#Here we are dealing with special situations in final results.
	if length == 0:
		ClusterResultWithUserID[PostID] = []
		ClusterResultWithUserID[PostID].append('0')

		continue


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
#	for w in range(0,200):
#		evals, Y, idx = cluster_points(L) 
	
#		membership = []

#		for i in range(0,len(idx)):
#			membership.append(str(idx[i]))
#			membership[i] = int(membership[i])

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
	print '--------------------' + PostID + '-----------------------'
	print 'modularity value is :' + str(checkmodularity)
	print 'OneGroup:'
	print OneGroup
	print 'AnotherGroup:'
	print AnotherGroup
#	print len(DegreeForOneGroup)
#	layout1 = gr.layout('fr')#fr forced directed algorithm
#	plot(gr,layout=layout1)
	break





