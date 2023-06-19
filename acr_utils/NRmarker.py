import operator
from ete3 import parser, Tree

def children(t, nodename):
    """
    find and rturn all leafs belpongin as chidlren to a node
    """
    try:
        nn = t.search_nodes(name=nodename)[0].get_leaf_names()
    except IndexError:
        print(f"Could not find any children for node {nodename}")
        nn = []
    return nn

def getSets(setinfo,ngfilter=3,nfilter=20): # from EukCC default: ngfilter3, nfilter=20
    sets=[]
    with open(setinfo,'r') as f:
        next(f)
        for i in f:
            n,ngenomes,node=i.rstrip('\n').split(',')
            if int(ngenomes)>=ngfilter and int(n)>=nfilter: # from EukCC default
                sets.append({'n':int(n),'ngenomes':int(ngenomes),'node':node})
    return sets

def find_N_marker(ori_tree,new_tree,setinfo):
    oriT = Tree(ori_tree, format=1)
    newT = Tree(new_tree,format=1)

    ori_names=oriT.get_tree_root().get_leaf_names()
    placement=set(newT.get_tree_root().get_leaf_names())-set(ori_names)

    setinfo=getSets(setinfo ,ngfilter=5,nfilter=50) #hoonje 20220912: ngfilter=5, nfilter=50
    tmp=[]
    for s in setinfo:
        cover=set(children(newT,s['node'])) & placement
        if cover:
            s["covering"]=list(placement)[0]
            tmp.append(s)
    setinfo=tmp
    if not setinfo:
        return None
    # top sorting
    # finding neareast marker set
    setinfo=sorted(setinfo, key=operator.itemgetter('ngenomes'), reverse=False) #HPA:reverse_True | LCA:reverse_False
    sisters=set()
    for s in newT.search_nodes(name=setinfo[0]['covering'])[0].get_sisters():
        sisters=sisters|set(s.get_leaf_names())

    #return a MAG, a nearest node, sisters
    return (setinfo[0]['covering'],setinfo[0]['node'],sisters)
