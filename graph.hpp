#pragma once
#include "dd.h"
#include <vector>
#include "common.h"
#include <atomic>
#include <iostream>
#include <thread>
#include <random>
#include <variant>
#include "wsq.hpp"

#ifdef GRAPHVIZ
#include "cgraph.h"
#include "gvc.h"
#endif

class Worker;

class Node {
    friend class Graph;
    public:

        struct IdentM{
            IdentM(const mEdge& m): work(
                std::bind([=](Worker* w){ return m;}, std::placeholders::_1)
                    ){}
            std::function<mEdge(Worker*)> work;
        };

        struct IdentV{
            IdentV(const vEdge& v): work(
                std::bind([=](Worker* w){ return v;}, std::placeholders::_1)
                    ){}
            std::function<vEdge(Worker*)> work;
        };


        struct AddMM{
            AddMM() = default;
            void bind(const mEdge& lhs, const mEdge& rhs){
                work = std::bind(mm_add, std::placeholders::_1, lhs, rhs);
            }
            std::function<mEdge(Worker*)> work; 
        };

        struct AddVV{
            AddVV() = default;
            void bind(const vEdge& lhs, const vEdge& rhs){
                work = std::bind(vv_add, std::placeholders::_1, lhs, rhs);        
            }
            std::function<vEdge(Worker*)> work; 
        };

        struct MulMM{
            MulMM() = default;
            void bind(const mEdge& lhs, const mEdge& rhs){
                work = std::bind(mm_multiply, std::placeholders::_1, lhs, rhs);
            }
            std::function<mEdge(Worker*)> work; 
        };

        struct MulMV{
            MulMV() = default;
            void bind(const mEdge& lhs, const vEdge& rhs){
                work = std::bind(mv_multiply, std::placeholders::_1, lhs, rhs);
            }
            std::function<vEdge(Worker*)> work; 
        };

        struct MulVV{
            MulVV() = default;
            void bind(const vEdge& lhs, const vEdge& rhs){
                work = std::bind(vv_multiply, std::placeholders::_1, lhs, rhs);
            }
            std::function<vEdge(Worker*)> work;

        };


        using Edge_t = std::variant<std::monostate,mEdge, vEdge>;
        constexpr static auto EDGE_PLACEHOLDER = get_index_v<std::monostate, Edge_t>;
        constexpr static auto MEDGE = get_index_v<mEdge, Edge_t>;
        constexpr static auto VEDGE = get_index_v<vEdge, Edge_t>;

        using task_t = std::variant<std::monostate, IdentM, IdentV, AddMM, AddVV, MulMM, MulMV, MulVV>;
        constexpr static auto TASK_PLACE_HOLDER = get_index_v<std::monostate, task_t>;
        constexpr static auto IDENTM = get_index_v<IdentM, task_t>;
        constexpr static auto IDENTV = get_index_v<IdentV, task_t>;
        constexpr static auto ADDMM = get_index_v<AddMM, task_t>;
        constexpr static auto ADDVV = get_index_v<AddVV, task_t>;
        constexpr static auto MULMM = get_index_v<MulMM, task_t>;
        constexpr static auto MULMV = get_index_v<MulMV, task_t>;
        constexpr static auto MULVV = get_index_v<MulVV, task_t>;

        
        Node(const mEdge& m): _required(0){
            task = IdentM(m);
        }

        Node(const vEdge& v): _required(0){
            task = IdentV(v);
        }

        Node(decltype(TASK_PLACE_HOLDER) t): _required(2), _dependents{2}{
            switch(t){
                case ADDMM:    task = AddMM();      break;
                case ADDVV:    task = AddVV();      break;
                case MULMM:    task = MulMM();      break;
                case MULMV:    task = MulMV();      break;
                case MULVV:    task = MulVV();      break;
                default:  {
                    std::cout<<"Unsupported task: "<<t<<std::endl;
                    exit(1);
                }
            }
        }

        void lhs(Node* v){
            assert(v->_dependents.size() == 2);
            _successors.push_back(v);
            v->_dependents[0] = this;
        }

        void rhs(Node* v){
            assert(v->_dependents.size() == 2);
            _successors.push_back(v);
            v->_dependents[1] = this;
        
        }

        size_t num_successors() const {
            return _successors.size();
        }

        size_t num_dependents() const {
            return _dependents.size();
        }

        bool ready_to_be_executed() const {
            return _joint.load(std::memory_order_acquire) == _required;
        }

        void prepare_to_be_executed() const {


                
            
        }
    
        void execute(Worker* w){
            std::size_t idx = task.index(); 
            switch(idx){
                case IDENTM:     update_result(w, std::get<IdentM>(task).work(w));  break;
                case IDENTV:     update_result(w, std::get<IdentV>(task).work(w));  break;
                case ADDMM:      update_result(w, std::get<AddMM>(task).work(w));   break;
                case ADDVV:      update_result(w, std::get<AddVV>(task).work(w));   break;
                case MULMM:      update_result(w, std::get<MulMM>(task).work(w));   break;
                case MULMV:      update_result(w, std::get<MulMV>(task).work(w));   break;
                case MULVV:      update_result(w, std::get<MulVV>(task).work(w));   break;
                default:         std::cout<<"Unsupported task"<<std::endl;       break;
            }
            
        }


        void update_result(Worker* w, const Edge_t& r){
            result = r;
            for(Node* s: _successors){
                s->_joint.fetch_add(1, std::memory_order_release);
                if(s->ready_to_be_executed()){
                    s->prepare_to_be_executed();
                }
            }
        }
        

    private:

        std::vector<Node*> _successors;
        std::vector<Node*> _dependents;

        task_t task;
        Edge_t result;

        const int _required;
        std::atomic<int> _joint{0}; 

#ifdef GRAPHVIZ
        Agnode_t* agn{nullptr};
        
#endif
        
};

class Graph {
    friend class Executor;
    public:

        template<typename... Args>
            Node* emplace(Args... args){
                Node* n = new Node(std::forward<Args>(args)...); 
                _nodes.push_back(n);
                return n;
            }

            void emplace(Node* n){
                _nodes.push_back(n);
            }

            Graph() = default;
            Graph(Graph&& other): _nodes(std::move(other._nodes)){}
            Graph(const Graph& other): _nodes(other._nodes){}
            Graph& operator=(Graph&& other){
                _nodes = std::move(other._nodes);
                return *this;
            }
            Graph& operator=(const Graph& other){
                _nodes = other._nodes;
                return *this;
            }

            std::size_t root_size() const {
                return _nodes.size();
            }

#ifdef GRAPHVIZ
        void set_node_label(Node* n){
            if(n->agn){
              switch(n->task.index()) {
                case Node::IDENTM:      agset(n->agn, "label", "IdentM");   break;
                case Node::IDENTV:      agset(n->agn, "label", "IdentV");   break;
                case Node::ADDMM:       agset(n->agn, "label", "AddMM");    break;
                case Node::ADDVV:       agset(n->agn, "label", "AddVV");    break;  
                case Node::MULMM:       agset(n->agn, "label", "MulMM");    break;
                case Node::MULMV:       agset(n->agn, "label", "MulMV");    break;
                case Node::MULVV:       agset(n->agn, "label", "MulVV");    break;
                default:                std::cout<<"unrecognized task"<<std::endl; exit(1); 
              }
            }
        }

        int dump() {
            Agraph_t* g;
            g = agopen("G", Agdirected, NULL);

            Agsym_t* sym = agattr(g, AGNODE, "label", "IdentM");
            for(Node* n: _nodes){
                Agnode_t* agn = agnode(g, NULL, true);
                if(agn == NULL){
                    std::cout<<"create graph node failed"<<std::endl;
                    exit(1);
                }
                n->agn = agn;
                set_node_label(n);

            }
            std::vector<Node*> current = _nodes;
            std::vector<Node*> next;
            while(!current.empty()){
                for(Node* n: current){
                    for(Node* head: n->_successors){
                        bool first_time = false;
                        if(head->agn == NULL){
                            head->agn = agnode(g, NULL, true);
                            set_node_label(head);
                            first_time = true;
                        }
                        Agedge_t* e = agedge(g, n->agn, head->agn, NULL, true);
                        if(first_time)next.push_back(head);
                    }

                }
                current.clear();
                current.swap(next);
            }
            agwrite(g, stdout);
            agclose(g);

            return 0;

        }
#else
        int dump(void* ofile){
            std::cout<<"Graphviz is not enabled"<<std::endl;
            return 1;
        }
#endif


    private:
        std::vector<Node*> _nodes; // only root nodes

};



class QuantumCircuit{
    public:
        QuantumCircuit(QubitCount q): _total_qubits(q){}

        template<typename... Args>
        void emplace_back(Args&&... args){
            _gates.emplace_back(makeGate(_total_qubits, std::forward<Args>(args)...));
        }


        void buildCircuit(){

            Graph g;
            std::vector<Node*> nodes;
            nodes.reserve(_gates.size());
            for(const mEdge& e: _gates){
                Node* n = new Node(e);
                nodes.push_back(n); 
                g.emplace(n);
            }

            mulmm_next_level(g, nodes, 0, nodes.size());
            _graph = std::move(g);
            return;

        }


        void dump_task_graph() {
            _graph.dump();
        }


    private:
        QubitCount _total_qubits;
        std::vector<mEdge> _gates;
        Graph _graph;
        vEdge _input;


        void mulmm_next_level(Graph& g,std::vector<Node*>& v, std::size_t start, std::size_t end){
            if(start == end - 1) return;
            std::size_t i;
            for(i = start; i < end - 1; i+=2){
                Node* n = new Node(Node::MULMM); 
                v[i]->lhs(n);
                v[i+1]->rhs(n);
                v.push_back(n);
            }

            if(i == end-1){
                v.push_back(v[i]);
            }

            mulmm_next_level(g,v, end, v.size());
        
        }

};

class Executor;

class Worker{
    friend class Executor;
    public:
    private:
        size_t _id;
        Executor* _executor;
        WorkStealingQueue<Node*> _wsq;
        std::thread* _thread;
        std::default_random_engine _rdgen { std::random_device{}() };
};

class Executor{
    public:
        Executor(int N): _nworkers(N), _workers{(std::size_t)N}{}
        void seed(const Graph& graph){
            for(Node* n: graph._nodes){
                _wsq.push(n);
            }
        }
        void spawn(); 
    private:
        void invoke(Worker*, Node*);
        void try_execute_self(Worker*);
        bool try_execute_else(Worker*);
        
        std::vector<Worker> _workers;
        WorkStealingQueue<Node*> _wsq;
        int _nworkers;
};
