// https://github.com/tskit-dev/tskit/blob/f27700635da6a5e20b885810f0f21a2f2203dcce/python/tests/__init__.py#L321
// implementation of schieber vishkin algorithm for finding the mrca of two nodes in a tree
#include <vector>
#include <iostream> // For std::cerr
#include <cstdlib>
#include "tskit.h"


#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }


// oriented forest: each tree will have its own oriented forest, which should be an iterable (vector in this case)
// that contains the parent of each node in the tree sequence for that tree

class MRCACalculator{

    public:
        MRCACalculator(const std::vector<int> &oriented_forest){
                std::vector<int> converted;
                converted.push_back(0);
                for(int i = 0; i < oriented_forest.size(); i++){
                    converted.push_back(oriented_forest[i]);
                }
                this->n = converted.size();
                preprocess(converted);
            }
        const int LAMBDA = 0;
        int get_mrca(int x, int y){
            return this->schieber_vishkin_mrca(x + 1, y + 1) - 1;
        };
    private:
        int n;
        std::vector<int> lambda;
        std::vector<int> pi;
        std::vector<int> tau;
        std::vector<int> beta;
        std::vector<int> alpha;
        
        

        void preprocess(std::vector<int> oriented_forest){
            std::vector<int> child(this->n, 0);
            std::vector<int> parent(this->n, 0);
            std::vector<int> sib(this->n, 0);
            bool notDone;
            for(int u = 0; u < this->n; u++){
                this->lambda.push_back(0);
                this->pi.push_back(0);
                this->tau.push_back(0);
                this->beta.push_back(0);
                this->alpha.push_back(0);
                int v = oriented_forest[u];
                sib[u] = child[v];
                child[v] = u;
                parent[u] = v; 
            }
            int p = child[this->LAMBDA];
            this->n= 0;
            this->lambda[0] = -1;
            while(p != this->LAMBDA){
                notDone = true;
                while(notDone){
                    this->n++;
                    this->pi[p] = n;
                    this->tau[n] = this->LAMBDA;
                    if(child[p] != this->LAMBDA){
                        p = child[p];
                    }
                    else{
                        notDone = false;
                    }
                }
                this->beta[p] = n;
                notDone = true;
                while(notDone){
                    this->tau[this->beta[p]] = parent[p];
                    if(sib[p] != this->LAMBDA){
                        p = sib[p];
                        notDone = false;
                    }
                    else{
                        p = parent[p];
                        if(p != this->LAMBDA){
                            int h = this->lambda[n & -this->pi[p]];
                            this->beta[p] = ((n >> h) | 1) << h; // huh? need to make sure signage is all good, since C++ and Python handle rightward bitwise shifts differently with negatives
                        }
                        else{
                            notDone = false;
                        }
                    }
                }
            }
            this->lambda[0] = this->lambda[n];
            this->pi[this->LAMBDA] = 0;
            this->beta[this->LAMBDA] = 0;
            this->alpha[this->LAMBDA] = 0;
            p = child[this->LAMBDA];
            while(p != this->LAMBDA){
                notDone = true;
                while(notDone){
                    int a = this->alpha[parent[p]] | (this->beta[p] & -this->beta[p]);
                    this->alpha[p] = a;
                    if(child[p] != this->LAMBDA){
                        p = child[p];
                    }
                    else{
                        notDone = false;
                    }
                }
                notDone = true;
                while(notDone){
                    if(sib[p] == this->LAMBDA){
                        p = sib[p];
                        notDone = false;
                    }
                    else{
                        p = parent[p];
                        notDone = (p != this->LAMBDA);
                    }
                }
            }
            
        }
        int schieber_vishkin_mrca(int x, int y){
            int h;
            if (this->beta[x] <= this->beta[y]){
                h = this->lambda[this->beta[y] & -this->beta[x]];
            }
            else{
                h = this->lambda[this->beta[x] & -this->beta[y]];
            }
            int k = this->alpha[x] & this->alpha[y] & -(1 << h);
            h = this->lambda[k & -k];
            int j = ((this->beta[x] >> h) | 1) << h;
            int xhat;
            int yhat;
            int ell;
            int z;
            if (j == this->beta[x]){
                xhat = x;
            }
            else{
                ell = this->lambda[this->alpha[x] & ((1 << h) - 1)];
                xhat = this->tau[((this->beta[x] >> ell) | 1) << ell];
            }
            if (j == this->beta[y]){
                yhat = y;
            }
            else{
                ell = this->lambda[this->alpha[y] & ((1 << h) - 1)];
                yhat = this->tau[((this->beta[y] >> ell) | 1) << ell];
            }
            if(this->pi[xhat] <= this->pi[yhat]){
                z = xhat;
            }
            else{
                z = yhat;
            }
            return z;
 
        }

};


int main(){
    char *ts_file = "../../data/gt_extraction_test.trees";

    tsk_treeseq_t ts;
    int ret = 0;
    ret = tsk_treeseq_load(&ts, ts_file, 0);
    check_tsk_error(ret);
    tsk_tree_t tree;
    check_tsk_error(ret);
    ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);
    ret = tsk_tree_last(&tree);
    check_tsk_error(ret);
    tsk_size_t n_nodes = tree.num_nodes;
    std::vector<int> parents;
    for (size_t i = 0; i < n_nodes; i++){
        tsk_id_t par;
        tsk_tree_get_parent(&tree, i, &par);
        int node = static_cast<int>(par);
        parents.push_back(node);
    }
    MRCACalculator calc(parents);
    int mrca_122_139 = calc.get_mrca(122, 139);
    std::cout << mrca_122_139 << std::endl;
    
}