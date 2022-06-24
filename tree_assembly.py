#构建测序碱基树，从而组装基因组
from Bio import SeqIO

class node():
    def __init__(self,count=1,this_base=None,level=0):
        self.count = count
        self.this_base = this_base
        self.level = level
        self.next_node()
        
    def next_node(self,A=None,T=None,C=None,G=None):
        self.A = A
        self.T = T
        self.C = C
        self.G = G

class tree(node):
    def __init__(self):
        self.A = None
        self.T = None
        self.C = None
        self.G = None
        self.this_base = None
        self.base_alias = {1:'A',2:'T',3:'C',4:'G',5:None,'A':1,'T':2,'C':3,'G':4}
    def make_tree(self,file):
        record = SeqIO(file,'fasta')
        for i in record:
            self.make_a_seq(str(i.seq))

    def make_a_seq(self,sequence):
        level = 1
        base_point = self
        for base in sequence:
            base_point = self.estimate_base(base,base_point,level)
            level += 1
            print(base_point.level)
        
    def estimate_base(self,base:str,base_point,level):
        base = base.upper()
        if base == 'A':
            base_point.A = self.update_point(base_point.A,'A',level)
            base_point = base_point.A
        elif base == 'T':
            base_point.T = self.update_point(base_point.T,'T',level)
            base_point = base_point.T
        elif base == 'C':
            base_point.C = self.update_point(base_point.C,'C',level)
            base_point = base_point.C
        elif base == 'G':
            base_point.G = self.update_point(base_point.G,'G',level)
            base_point = base_point.G
        else :
            raise InputError('the nucleobase is not ture base in the sequence')
        return base_point

    def update_point(self,next_point,base,level):
        if next_point:
            next_point.count += 1
        else:
            next_point = node(count=1,this_base=base,level=level)
        return next_point
    
    def next_node(self,a_node,next):
        #序列ATCG
        if next == 'A' :
            if a_node.A:
                return a_node.A
            else:
                next ='T'
        if next == 'T' :
            if a_node.T:
                return a_node.T
            else:
                next ='C'
        if next == 'C':
            if a_node.C:
                return a_node.C
            else :
                next ='G'
        if next == 'G' :
            if a_node.G:
                return a_node.G
            else:
                return False
        else:
            False 

    def end_dispose(self,end_node):
        end_base = end_node.this_base
        alias = self.base_alias[end_base]
        next_base = self.base_alias[alias+1]
        return next_base


    def before_node(self,node_list): 
        next_nodes = False
        while not len(node_list) == 1:
            end_node = node_list[-1]
            node_list.pop()
            next_base = self.end_dispose(end_node)
            base = next_base
            node_point = node_list[-1]
            next_nodes = self.next_node(node_point,base)
            if next_nodes:
                node_list.append(next_nodes)
                return node_list
        return node_list

    def get_seqs(self):
        node_list = [self]
        base = 'A'
        while node_list:
            node_point = node_list[-1]
            next_nodes = self.next_node(node_point,base)
            if next_nodes:
                node_list.append(next_nodes)
                base = 'A'
                continue
            else:

                yield node_list[1:]
                node_list =  self.before_node(node_list)
                if node_list[-1] == self:
                    break
                
        return 1

def test():
    seqs = 'Agctagcatgtcagtagcat'
    seqs2 = 'gctagctagcatgcatg'
    seqs3 = 'gctagctagcgtagct'
    a = tree()
    a.make_a_seq(seqs)
    a.make_a_seq(seqs2)
    a.make_a_seq(seqs3)
    for i in a.get_seqs():
        print([j.this_base for j in i])

if __name__ == '__main__':
    test()