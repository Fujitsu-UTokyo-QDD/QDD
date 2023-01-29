#pragma once
#include "dd.h"
#include <vector>
#include "common.h"
#include <atomic>
#include <iostream>
#include <thread>
#include <random>
#include <variant>
#include <semaphore>
#include "cache.hpp"
#include "wsq.hpp"
#include <cmath>
#include <semaphore>

#ifdef GRAPHVIZ
#include "cgraph.h"
#include "gvc.h"
#endif

#pragma GCC diagnostic ignored "-Wwrite-strings"

inline unsigned findPreviousPowerOf2(unsigned n)
{
    // drop all set bits from `n` except its last set bit
    return 1U << (int)log2(n);
}

class Node;
class Executor;

struct Worker{
        Worker(QubitCount q): _addCache(q), _mulCache(q){}
        

        void push_to_next_round(Node* n) {next_round.push_back(n);}
        void push_to_this_round(Node* n) {this_round.push_back(n);}

        size_t _id;
        Executor* _executor;
        std::thread* _thread{nullptr};
        std::default_random_engine _rdgen { std::random_device{}() };
        std::uniform_int_distribution<int> _dist;
        WorkStealingQueue<Node*> _wsq;

        AddCache _addCache;
        MulCache _mulCache;


        std::vector<Node*> this_round;
        std::vector<Node*> next_round;

        void execute();

        void collect(int N);
};

class Node {
    friend class Graph;
    friend class QuantumCircuit;
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
        struct ReduceMHub{
            ReduceMHub() = default;
           void bind (const mEdge& m){
                work = std::bind([=](Worker* w){ return m;}, std::placeholders::_1);
            }
            std::function<mEdge(Worker*)> work;
        };

        struct ReduceVHub{
            ReduceVHub() = default;
            void bind(const vEdge& v){ 
                work=std::bind([=](Worker* w){ return v;}, std::placeholders::_1);
            }
            std::function<vEdge(Worker*)> work;
        };
        struct ReduceMRep{
            ReduceMRep() = default;
           void bind (const mEdge& m){
                work = std::bind([=](Worker* w){ return m;}, std::placeholders::_1);
            }
            std::function<mEdge(Worker*)> work;
        };

        struct ReduceVRep{
            ReduceVRep() = default;
            void bind(const vEdge& v){ 
                work=std::bind([=](Worker* w){ return v;}, std::placeholders::_1);
            }
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

        using task_t = std::variant<std::monostate, IdentM, IdentV, AddMM, AddVV, MulMM, MulMV, MulVV, ReduceMHub, ReduceVHub, ReduceMRep, ReduceVRep>;
        constexpr static auto TASK_PLACE_HOLDER = get_index_v<std::monostate, task_t>;
        constexpr static auto IDENTM = get_index_v<IdentM, task_t>;
        constexpr static auto IDENTV = get_index_v<IdentV, task_t>;
        constexpr static auto ADDMM = get_index_v<AddMM, task_t>;
        constexpr static auto ADDVV = get_index_v<AddVV, task_t>;
        constexpr static auto MULMM = get_index_v<MulMM, task_t>;
        constexpr static auto MULMV = get_index_v<MulMV, task_t>;
        constexpr static auto MULVV = get_index_v<MulVV, task_t>;
        constexpr static auto REDUCEMHUB = get_index_v<ReduceMHub, task_t>;
        constexpr static auto REDUCEVHUB = get_index_v<ReduceVHub, task_t>;
        constexpr static auto REDUCEMREP = get_index_v<ReduceMRep, task_t>;
        constexpr static auto REDUCEVREP = get_index_v<ReduceVRep, task_t>;

        
        Node(const mEdge& m): _required(0){
            _task = IdentM(m);
        }

        Node(const vEdge& v): _required(0){
            _task = IdentV(v);
        }

        Node(Node* dep): _required(1){
            if(std::holds_alternative<Node::MulMV>(dep->_task)){
                _task = ReduceVHub();
                _dependents.push_back(dep);
                dep->_successors.push_back(this);
            }else if(std::holds_alternative<Node::MulMM>(dep->_task)){
                _task = ReduceMHub();
                _dependents.push_back(dep);
                dep->_successors.push_back(this);
            }else if (std::holds_alternative<Node::ReduceMHub>(dep->_task)){
                _task = ReduceMRep();
                _dependents.push_back(dep);
                dep->_successors.push_back(this);
            }else if (std::holds_alternative<Node::ReduceVHub>(dep->_task)){
                _task = ReduceVRep();
                _dependents.push_back(dep);
                dep->_successors.push_back(this);
            }else if (std::holds_alternative<Node::IdentM>(dep->_task)){
                _task = ReduceMHub();
                _dependents.push_back(dep);
                dep->_successors.push_back(this);
            }else if (std::holds_alternative<Node::IdentV>(dep->_task)){
                _task = ReduceVHub();
                _dependents.push_back(dep);
                dep->_successors.push_back(this);
            }else{
                std::cout<<"Unsupported reduce task"<<std::endl;
                exit(1);
            }
        }



        Node(decltype(TASK_PLACE_HOLDER) t): _required(2), _dependents{2}{
            switch(t){
                case ADDMM:    _task = AddMM();      break;
                case ADDVV:    _task = AddVV();      break;
                case MULMM:    _task = MulMM();      break;
                case MULMV:    _task = MulMV();      break;
                case MULVV:    _task = MulVV();      break;
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

        void precede(Node* v){
            _successors.push_back(v);
            v->_dependents.push_back(this);
            v->_required++;
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
        

        void prepare_to_be_executed() {
            std::size_t idx = _task.index(); 
            switch(idx){
                case ADDMM:      std::get<AddMM>(_task).bind(std::get<MEDGE>(_dependents[0]->_result), std::get<MEDGE>(_dependents[1]->_result));   break;
                case ADDVV:      std::get<AddVV>(_task).bind(std::get<VEDGE>(_dependents[0]->_result), std::get<VEDGE>(_dependents[1]->_result));   break;
                case MULMM:      std::get<MulMM>(_task).bind(std::get<MEDGE>(_dependents[0]->_result), std::get<MEDGE>(_dependents[1]->_result));   break;
                case MULMV:      std::get<MulMV>(_task).bind(std::get<MEDGE>(_dependents[0]->_result), std::get<VEDGE>(_dependents[1]->_result));   break;
                case MULVV:      std::get<MulVV>(_task).bind(std::get<VEDGE>(_dependents[0]->_result), std::get<VEDGE>(_dependents[1]->_result));   break;
                case REDUCEMHUB:    std::get<ReduceMHub>(_task).bind(std::get<MEDGE>(_dependents[0]->_result));                                     break;
                case REDUCEVHUB:    std::get<ReduceVHub>(_task).bind(std::get<VEDGE>(_dependents[0]->_result));                                     break;
                case REDUCEMREP:    std::get<ReduceMRep>(_task).bind(std::get<MEDGE>(_dependents[0]->_result));                                     break;
                case REDUCEVREP:    std::get<ReduceVRep>(_task).bind(std::get<VEDGE>(_dependents[0]->_result));                                     break;
                case IDENTM:     
                case IDENTV:     break;
                default:         std::cout<<"Unsupported task"<<std::endl;       break;
            }


                
            
        }
    
        void execute(Worker* w){
            std::size_t idx =_task.index(); 
            switch(idx){
                case IDENTM:     update_result(w, std::get<IdentM>(_task).work(w));  break;
                case IDENTV:     update_result(w, std::get<IdentV>(_task).work(w));  break;
                case ADDMM:      update_result(w, std::get<AddMM>(_task).work(w));   break;
                case ADDVV:      update_result(w, std::get<AddVV>(_task).work(w));   break;
                case MULMM:      update_result(w, std::get<MulMM>(_task).work(w));   break;
                case MULMV:      update_result(w, std::get<MulMV>(_task).work(w));   break;
                case MULVV:      update_result(w, std::get<MulVV>(_task).work(w));   break;
                case REDUCEMHUB:    update_result(w, std::get<ReduceMHub>(_task).work(w)); break;
                case REDUCEVHUB:    update_result(w, std::get<ReduceVHub>(_task).work(w)); break;
                case REDUCEMREP:    update_result(w, std::get<ReduceMRep>(_task).work(w)); break;
                case REDUCEVREP:    update_result(w, std::get<ReduceVRep>(_task).work(w)); break;
                default:         std::cout<<"Unsupported task"<<std::endl;       break;
            }
            
        }


        void update_result(Worker* w, const Edge_t& r){
            assert(_result.index() == EDGE_PLACEHOLDER);
            _result = r;
            for(Node* s: _successors){
                if(s->_joint.fetch_add(1, std::memory_order_release) == (s->_required-1)){
                    s->prepare_to_be_executed();
                    w->push_to_next_round(s);
                }
            }
            if(this->_sem != nullptr) this->_sem->release();
        }
        

    private:

        std::vector<Node*> _successors;
        std::vector<Node*> _dependents;

        task_t _task;
        Edge_t _result;

        int _required;
        std::atomic<int> _joint{0}; 
        std::binary_semaphore* _sem{nullptr};

#ifdef GRAPHVIZ
        Agnode_t* agn{nullptr};
        
#endif
        
};

class Graph {
    friend class Executor;
    friend class QuantumCircuit;
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
              switch(n->_task.index()) {
                case Node::IDENTM:      agset(n->agn, "label", "IdentM");   break;
                case Node::IDENTV:      agset(n->agn, "label", "IdentV");   break;
                case Node::ADDMM:       agset(n->agn, "label", "AddMM");    break;
                case Node::ADDVV:       agset(n->agn, "label", "AddVV");    break;  
                case Node::MULMM:       agset(n->agn, "label", "MulMM");    break;
                case Node::MULMV:       agset(n->agn, "label", "MulMV");    break;
                case Node::MULVV:       agset(n->agn, "label", "MulVV");    break;
                case Node::REDUCEMHUB:     agset(n->agn, "label", "ReduceMHub");    break;
                case Node::REDUCEVHUB:     agset(n->agn, "label", "ReduceVHub");    break;
                case Node::REDUCEMREP:     agset(n->agn, "label", "ReduceMRep");    break;
                case Node::REDUCEVREP:     agset(n->agn, "label", "ReduceVRep");    break;
                default:                std::cout<<"unrecognized task"<<std::endl; exit(1); 
              }
            }
        }

        int dump() {
            Agraph_t* g;
            g = agopen("G", Agdirected, NULL);

            Agsym_t* node_sym = agattr(g, AGNODE, "label", "IdentM");
            Agsym_t* edge_sym = agattr(g, AGEDGE, "label", "L");
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
                    for(int i = 0; i < n->_successors.size(); i++){
                        Node* head = n->_successors[i];
                        bool first_time = false;
                        if(head->agn == NULL){
                            head->agn = agnode(g, NULL, true);
                            set_node_label(head);
                            first_time = true;
                        }
                        Agedge_t* e = agedge(g, n->agn, head->agn, NULL, true);
                        if(std::holds_alternative<Node::ReduceMHub>(head->_task) || std::holds_alternative<Node::ReduceVHub>(head->_task) || std::holds_alternative<Node::ReduceMRep>(head->_task)||std::holds_alternative<Node::ReduceVRep>(head->_task)||std::holds_alternative<Node::IdentM>(head->_task)||std::holds_alternative<Node::IdentV>(head->_task)){
                                agset(e, "label", "");
                        }else if(n == head->_dependents[1]){
                            agset(e, "label", "R");
                        }

                        
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
        int dump(){
            std::cout<<"Graphviz is not enabled"<<std::endl;
            return 1;
        }
#endif


    private:
        std::vector<Node*> _nodes; // only root nodes

};


class Executor{
    friend class Worker;
    public:
        Executor(int N, bool* s, QubitCount q): _nworkers(N), _stop(s), _sem(0){
            for(int i = 0; i < N; i++) _workers.emplace_back(new Worker(q));
            _total_queue.resize(N);
        }

        ~Executor(){
            for(Worker* w: _workers){
                if(w->_thread != nullptr){
                    w->_thread->join();
                }
            }

            duration_micro longest = steal_timer.combine([](const duration_micro& t1, const duration_micro& t2){
                return (t1 > t2)? t1: t2;
            });
            duration_micro sum = steal_timer.combine([](const duration_micro& t1, const duration_micro& t2){
                return t1 + t2 ;
            });

            //std::cout<<"longest time spent in wait: "<<longest.count()<<" ms"<<std::endl;
            //std::cout<<"avg time spent in wait: "<<sum.count()/_nworkers<<" ms"<<std::endl;
        }

        void seed(const Graph& graph){
            int ntasks = graph._nodes.size() / _nworkers;
            std::cout<<"ntasks before:"<<ntasks<<std::endl;
            ntasks = findPreviousPowerOf2(ntasks);
            std::cout<<"ntasks after:"<<ntasks<<std::endl;
            int w = 0;
            for(;w < _nworkers; w++){
                for(int t = 0; t < ntasks; t++){
                    _workers[w]->push_to_this_round(graph._nodes[w*ntasks+t]);
                }
            }

            for(int t = w*ntasks; t < graph._nodes.size(); t++ )
                _workers[w-1]->push_to_this_round(graph._nodes[t]);
        }


        void spawn(); 
    private:
        const std::size_t MAX_STEALS{100};
        void try_execute_self(Worker*);
        void try_execute_else(Worker*);
        
        std::vector<Worker*> _workers;
        int _nworkers;
        bool* _stop;

        std::vector<Node*> _total_queue;
        std::counting_semaphore<64> _sem;

};

class QuantumCircuit{
    public:
        QuantumCircuit(QubitCount q, int nworkers, int reduce): _total_qubits(q), _stop(false), REDUCE_THRESHOLD(reduce), _executor(nworkers, &_stop, q){}
        QuantumCircuit(QubitCount q, int nworkers, int reduce, const vEdge& v): _total_qubits(q), _stop(false),REDUCE_THRESHOLD(reduce), _executor(nworkers, &_stop, q), _input(v){}

        template<typename... Args>
        void emplace_back(Args&&... args){
            _gates.emplace_back(makeGate(_total_qubits, std::forward<Args>(args)...));
        }

        void emplace_gate(const mEdge& e){
            _gates.emplace_back(e);
        }

        void emplace_add(){
            add_pos.emplace_back(_gates.size());
        }

        /*
        void buildCircuit(){

            std::vector<Node*> nodes;
            
            if(_input.n != nullptr){
                Node* n = new Node(_input);
                nodes.push_back(n);
            }
            
            std::size_t i = 0;
            for(; i < _gates.size()&& i < REDUCE_THRESHOLD; i++){
                Node* n = new Node(_gates[i]);
                nodes.push_back(n);
            }

            for(Node* n: nodes) _graph.emplace(n);
            
            Node* root = reduce_nodes_with_mul(nodes);  
            root = mul_next_level_reduce(root, _gates, _gates.size(), i, i + REDUCE_THRESHOLD);
            _executor.seed(_graph);
            assert(root != nullptr);

            _output = root;
            root->_sem = new std::binary_semaphore(0);
            return;

        }
        */
        void buildCircuit(){

            Graph g;
            std::vector<Node*> nodes;
            nodes.reserve(_gates.size() + 1);
            
            if(_input.n != nullptr){
                Node* n = new Node(_input);
                nodes.push_back(n);
                g.emplace(n);
            }

            for(const mEdge& e: _gates){
                Node* n = new Node(e);
                nodes.push_back(n); 
                g.emplace(n);
            }

            mul_next_level(g, nodes, 0, nodes.size());
            _graph = std::move(g);
            _executor.seed(_graph);
            return;

        }

        void setInput(const vEdge& v){ _input = v; }



        QuantumCircuit& wait(){
            _executor.spawn();
            assert(_output->_sem != nullptr);
            _output->_sem->acquire();
           _stop = true; 
           return *this;
        }

        mEdge matrixResult() {
            return std::get<Node::MEDGE>(_output->_result);
        }

        vEdge vectorResult(){
            return std::get<Node::VEDGE>(_output->_result);
        }


        void dump_task_graph() {
            _graph.dump();
        }

        QubitCount getQubits() const{
            return _total_qubits; 
        }

        


    private:

        std::size_t REDUCE_THRESHOLD{100};

        QubitCount _total_qubits;
        std::vector<mEdge> _gates;
        std::vector<int> add_pos; // for example, if add_pos == 3, this means we have 0,1, 2, in total 3 gates multiplied together, followed by an addition


        Graph _graph;
        bool _stop;
        vEdge _input;
        Node* _output;

        Executor _executor;


        void mul_next_level(Graph& g,std::vector<Node*>& v, std::size_t start, std::size_t end){
            if(start == end - 1) {
                v[start]->_sem = new std::binary_semaphore(0);
                _output = v[start];
                return;
            }
            std::size_t i;
            for(i = start; i < end - 1; i+=2){
                Node* n;
                if(std::holds_alternative<Node::IdentV>(v[i]->_task) || std::holds_alternative<Node::MulMV>(v[i]->_task)){
                    n = new Node(Node::MULMV);
                }else{
                     n = new Node(Node::MULMM); 
                }
                v[i]->rhs(n);
                v[i+1]->lhs(n);
                v.push_back(n);
            }

            if(i == end-1){
                v.push_back(v[i]);
            }

            mul_next_level(g,v, end, v.size());
        
        }
        
        void brp(){}

        Node* reduce_nodes_with_mul(std::vector<Node*> root /*we need a copy*/){
            std::vector<Node*> next_root;

            
            while( root.size() > 1 ){
                for(auto i = 0; i < root.size()/2; i++){
                    Node* n;
                    if(std::holds_alternative<Node::IdentV>(root[i]->_task) || std::holds_alternative<Node::MulMV>(root[i]->_task) || std::holds_alternative<Node::ReduceVRep>(root[i]->_task)){
                        n = new Node(Node::MULMV);
                    }else{
                        n = new Node(Node::MULMM); 
                    }
                    root[i*2]->rhs(n);
                    root[i*2+1]->lhs(n);
                    next_root.emplace_back(n);
                }

                if(root.size() % 2 == 1){
                    next_root.emplace_back(root.back());
                }

                root.clear();
                root.swap(next_root);
            }
            
            if(root.size() == 1)
                return new Node(root.back());
            else 
            {
                brp();
                return nullptr;
            }
            
            
        }


        Node* mul_next_level_reduce(Node* root, const std::vector<mEdge>& gates, std::size_t gates_end, std::size_t start, std::size_t end){

            assert(std::holds_alternative<Node::ReduceVHub>(root->_task) || std::holds_alternative<Node::ReduceMHub>(root->_task));
            
            std::vector<Node*> level;

            if(start < gates_end){
                Node* n = new Node(root);
                level.push_back(n);
            }

            if(start >= gates_end){
                return root;
            }

            std::size_t i = start;

            for(; i < end && i < gates_end; i++){
                Node* n = new Node(gates[i]);
                if(root != nullptr)
                    root->precede(n);
                level.push_back(n);
            }

            root = reduce_nodes_with_mul(level);


            return mul_next_level_reduce(root, gates, gates_end, end, end + REDUCE_THRESHOLD);

        
        }

};


