//
// Created by wouter on 29-10-19.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#define NBS_INIT_SIZE 32
#define TREE_INIT_SIZE 32
#define MAX_DIAMETER 1000


struct simple_doubly_linked_list_node {
    struct simple_doubly_linked_list_node* next;
    struct simple_doubly_linked_list_node* prev;
    unsigned int data;
};

struct simple_doubly_linked_list {
    struct simple_doubly_linked_list_node* first;
    struct simple_doubly_linked_list_node* last;
};

void init_simple_linked_list (struct simple_doubly_linked_list* ll) {
    ll->first = NULL;
    ll->last = NULL;
};

void simple_append(struct simple_doubly_linked_list* sll, unsigned int data) {
    struct simple_doubly_linked_list_node *new_node = malloc(sizeof(struct simple_doubly_linked_list_node));
    new_node->data = data;
    new_node->next = NULL;
    if (sll->first == NULL) {
        sll->first = new_node;
        sll->last = new_node;
        new_node->prev = NULL;
    } else {
        sll->last->next = new_node;
        new_node->prev = sll->last;
        sll->last = new_node;
    }
}

unsigned int simple_pop(struct simple_doubly_linked_list* ll) {
    if (ll->first == NULL) {
        return 0;
    }
    struct simple_doubly_linked_list_node *front = ll->first;
    unsigned int ret = front->data;
    ll->first = front->next;
    if (ll->first != NULL) {
        ll->first->prev = NULL;
    } else {
        ll->last = NULL;
    }
    free(front);
    return ret;
}

unsigned int simple_pop_right(struct simple_doubly_linked_list* ll) {
    if (ll->first == NULL) {
        return 0;
    }
    struct simple_doubly_linked_list_node *back = ll->last;
    unsigned int ret = back->data;
    ll->last = back->prev;
    if (ll->last != NULL) {
        ll->last->next = NULL;
    } else {
        ll->first = NULL;
    }
    free(back);
    return ret;
}

bool simple_is_empty(struct simple_doubly_linked_list* ll) {
    return ll->first == NULL;
}


struct linked_list_node {
    struct linked_list_node* next;
    void* data;
};

struct linked_list {
    struct linked_list_node* first;
    struct linked_list_node* last;
};

void init_linked_list (struct linked_list* ll) {
    ll->first = NULL;
    ll->last = NULL;
};

void append(struct linked_list* ll, void* data) {
    struct linked_list_node *new_node = malloc(sizeof(struct linked_list_node));
    new_node->data = data;
    new_node->next = NULL;
    if (ll->first == NULL) {
        ll->first = new_node;
        ll->last = new_node;
    } else {
        ll->last->next = new_node;
        ll->last = new_node;
    }
}

void* pop(struct linked_list* ll) {
    if (ll->first == NULL) {
        return NULL;
    }
    struct linked_list_node *front = ll->first;
    void* ret = front->data;
    ll->first = front->next;
    free(front);
    return ret;
}

bool is_empty(struct linked_list* ll) {
    return ll->first == NULL;
}


struct nd_tuple {
    unsigned int node;
    unsigned int dist;
};

struct node {
    unsigned int nbs_size;
    unsigned int max_nbs_size;
    unsigned int* nbs_list;
    unsigned int parent;
    unsigned int tree_depth;
    unsigned int max_tree_depth;
    unsigned int* tree_struct;
    unsigned int num_children;
};

struct edge {
    unsigned int head;
    unsigned int tail;
};

struct es_graph {
    unsigned int vsize;
    unsigned int esize;
    double d;
    struct node* v;
    struct edge* e;
    unsigned int* distribution;
    unsigned int* component_sizes;
    unsigned int* component_tree;
    unsigned long total_distance;
    unsigned long total_distance_div;
    unsigned long* typical_distance;
    unsigned int diameter;
};

//double rand_double() {
//    return (double) random() / (double) (RAND_MAX);
//}

double rand_double() {
    return drand48();
}


/* Returns an integer in the range [0, n).
 *
 * Uses rand(), and so is affected-by/affects the same seed.
 */
long rand_int(unsigned int n) {
    if ((n - 1) == RAND_MAX) {
        return lrand48();
    } else {
        // Supporting larger values for n would requires an even more
        // elaborate implementation that combines multiple calls to rand()
        assert (n <= RAND_MAX);

        // Chop off all of the values that would cause skew...
        unsigned int end = RAND_MAX / n; // truncate skew
        assert (end > 0);
        end *= n;

        // ... and ignore results from rand() that fall above that limit.
        // (Worst case the loop condition should succeed 50% of the time,
        // so we can expect to bail out of this loop pretty quickly.)
        long r;
        while ((r = lrand48()) >= end);

        return r % n;
    }
}

//long rand_int(long max) {
//    return random() % max;
//}

void evolve(struct es_graph *g);

void build_nodes(struct es_graph *g);

void collapse_neighbors(struct es_graph *g);

double p();

unsigned long total_distance2(struct es_graph *g);

struct node* init_node(struct node *n) {
    n->max_nbs_size = NBS_INIT_SIZE;
    n->max_tree_depth = TREE_INIT_SIZE;
    n->tree_depth = 1;
    n->nbs_size = 0;
    n->tree_struct = calloc(n->max_tree_depth, sizeof(int));
    n->tree_struct[0] = 1;
    n->nbs_list = malloc(sizeof(int)* n->max_nbs_size);
    n->num_children = 0;
    return n;
}

struct es_graph *generate_edge_step(double d, int t) {
    struct es_graph *g = malloc(sizeof(struct es_graph));
    g->vsize = 1;
    g->esize = 1;
    g->d = d;
    g->v = (struct node*)calloc(t, sizeof(struct node));
    g->e = (struct edge*)malloc(sizeof(struct edge)*t);
    g->distribution = (unsigned int*)malloc(sizeof(int)*t*2);
    g->component_sizes = (unsigned int*)malloc(sizeof(int)*t);
    g->component_tree = (unsigned int*)malloc(sizeof(int)*t);
    g->total_distance = 0;
    g->total_distance_div = 0;
    g->diameter = 0;
    g->typical_distance = (unsigned long*)calloc(MAX_DIAMETER, sizeof(long));

    g->component_sizes[0] = 1;
    g->component_tree[0] = 0;
    g->e[0].head = 0;
    g->e[0].tail = 0;
    g->distribution[0] = 0;
    g->distribution[1] = 0;
    init_node(&g->v[0]);

    for (unsigned int i = 1; i < t; i++) {
        evolve(g);
        printf("%u: %lu, %f\n", i, g->total_distance, (double)(2*g->total_distance)/(double)g->total_distance_div);
    }
//    build_nodes(g);
//    collapse_neighbors(g);
    return g;
}

void collapse_neighbors(struct es_graph *g) {
    bool* encounters = (bool*)calloc(g->vsize, sizeof(bool));
    memset(encounters, false, g->vsize*sizeof(bool));
    for (int i = 0; i < g->vsize; i++) {
        for (int j = 0; j < g->v[i].nbs_size; j++) {
            encounters[g->v[i].nbs_list[j]] = true;
        }
        encounters[i] = false;
        unsigned int tmp = 0;
        for (int j = 0; j < g->v[i].nbs_size; j++) {
            if (encounters[g->v[i].nbs_list[j]]) {
                g->v[i].nbs_list[tmp] = g->v[i].nbs_list[j];
                tmp++;
                encounters[g->v[i].nbs_list[j]] = false;
            }
        }
        g->v[i].nbs_size = tmp;
    }
    free(encounters);
}

void build_nodes(struct es_graph *g) {
    for (unsigned int i = 0; i < g->esize; i++) {
        g->v[g->e[i].head].nbs_size++;
        g->v[g->e[i].tail].nbs_size++;
    }
    for (unsigned int i = 0; i < g->vsize; i++) {
        struct node *n = &(g->v[i]);
        n->nbs_list = (unsigned int*)malloc(sizeof(unsigned int*)*n->nbs_size);
        n->nbs_size = 0;
    }
    for (unsigned int i = 0; i < g->esize; i++) {
        unsigned int head = g->e[i].head;
        unsigned int tail = g->e[i].tail;
        g->v[head].nbs_list[g->v[head].nbs_size] = tail;
        g->v[head].nbs_size++;
        g->v[tail].nbs_list[g->v[tail].nbs_size] = head;
        g->v[tail].nbs_size++;
    }
}

unsigned int find_component_root(struct es_graph *g, unsigned int n) {
    while (g->component_tree[n] != n) {
        n = g->component_tree[n];
    }
    return n;
}


void resize_tree(struct node *n) {
    unsigned int old_size = n->max_tree_depth;
    n->max_tree_depth *= 2;
    n->tree_struct = realloc(n->tree_struct, sizeof(int)*n->max_tree_depth);
    memset(n->tree_struct+sizeof(int)*old_size, 0, sizeof(int)*old_size);
}

void resize_nbs(struct node *n) {
    n->max_nbs_size *= 2;
    n->nbs_list = realloc(n->nbs_list, sizeof(int)*n->max_nbs_size);
}


void update_distance(struct es_graph *g) {
    bool* found = malloc(sizeof(bool)*g->vsize);
    struct nd_tuple* q = malloc(sizeof(struct nd_tuple)*g->vsize);
    unsigned int qhead = 0;
    unsigned int qtail = 0;

    unsigned long total = 0;

    for (unsigned int i = 0; i < g->vsize; i++) {
        if (g->v[i].parent != i) continue;
        memset(found, false, sizeof(bool)*g->vsize);
        found[i] = true;
        qhead = 0;
        qtail = 0;
        q[qtail].node = i;
        q[qtail].dist = 0;
        qtail++;


        while (qhead != qtail) {
            struct nd_tuple cur = q[qhead];
            qhead++;
            total += cur.dist * (g->v[cur.node].num_children+1) * (g->v[i].num_children+1);
            for (unsigned int j = 0; j < g->v[cur.node].nbs_size; j++) {
                unsigned int nb_id = g->v[cur.node].nbs_list[j];
                if (!found[nb_id] && g->v[nb_id].parent == nb_id) {
                    q[qtail].node = nb_id;
                    q[qtail].dist = cur.dist + 1;
                    qtail++;
                    found[nb_id] = true;
                }
            }
        }

    }
    total >>= 1u;

    for (unsigned int i = 0; i < g->vsize; i++) {
        if (g->v[i].parent == i) continue;
        total += (g->v[i].num_children+1) * (g->component_sizes[find_component_root(g, i)] - (g->v[i].num_children+1));
    }
    free(q);
    free(found);
    g->total_distance = total;
}


void edge_step(struct es_graph *g) {
    unsigned int r1, r2;
    if (rand_double() * (2*g->esize + 2*g->d) >= 2*g->d) {
        r1 = g->distribution[rand_int(2*g->esize)];
    } else {
        r1 = (unsigned int)rand_int(g->vsize);
    }
    if (rand_double() * (2*g->esize + 2*g->d) >= 2*g->d) {
        r2 = g->distribution[rand_int(2*g->esize)];
    } else {
        r2 = (unsigned int)rand_int(g->vsize);
    }
    struct edge *e = &g->e[g->esize];
    e->head = r1;
    e->tail = r2;
    g->distribution[2*g->esize] = r1;
    g->distribution[2*g->esize+1] = r2;

    bool exists = (r1 == r2);
    if (!exists) {
        if (g->v[r1].nbs_size < g->v[r2].nbs_size) {
            for (int i = 0; i < g->v[r1].nbs_size; i++) {
                if (g->v[r1].nbs_list[i] == r2) {
                    exists = true;
                    break;
                }
            }
        } else {
            for (int i = 0; i < g->v[r2].nbs_size; i++) {
                if (g->v[r2].nbs_list[i] == r1) {
                    exists = true;
                    break;
                }
            }
        }
    }
    if (!exists) {
        struct simple_doubly_linked_list sdll;
        init_simple_linked_list(&sdll);
        simple_append(&sdll, r1);
        unsigned int tmp = r1;
        while (g->v[tmp].parent != tmp) {
            unsigned int parent = g->v[tmp].parent;
            simple_append(&sdll, parent);
            g->v[tmp].parent = tmp;
            tmp = parent;
        }
        unsigned int prev = simple_pop_right(&sdll);
        while (!simple_is_empty(&sdll)) {
            unsigned int next = simple_pop_right(&sdll);
            for (unsigned int i = 0; i < g->v[next].tree_depth; i++) {
                g->v[prev].tree_struct[i + 1] -= g->v[next].tree_struct[i];
                g->v[prev].num_children -= g->v[next].tree_struct[i];
            }
            while (g->v[prev].tree_struct[g->v[prev].tree_depth-1] == 0) {
                g->v[prev].tree_depth--;
            }
            prev = next;
        }
        simple_append(&sdll, r2);
        tmp = r2;
        while (g->v[tmp].parent != tmp) {
            unsigned int parent = g->v[tmp].parent;
            simple_append(&sdll, parent);
            g->v[tmp].parent = tmp;
            tmp = parent;
        }
        prev = simple_pop_right(&sdll);
        while (!simple_is_empty(&sdll)) {
            unsigned int next = simple_pop_right(&sdll);
            for (unsigned int i = 0; i < g->v[next].tree_depth; i++) {
                g->v[prev].tree_struct[i + 1] -= g->v[next].tree_struct[i];
                g->v[prev].num_children -= g->v[next].tree_struct[i];
            }
            prev = next;
        }

        g->v[r1].nbs_size += 1;
        if (g->v[r1].nbs_size > g->v[r1].max_nbs_size) {
            resize_nbs(&g->v[r1]);
        }
        g->v[r1].nbs_list[g->v[r1].nbs_size-1] = r2;

        g->v[r2].nbs_size += 1;
        if (g->v[r2].nbs_size > g->v[r2].max_nbs_size) {
            resize_nbs(&g->v[r2]);
        }
        g->v[r2].nbs_list[g->v[r2].nbs_size-1] = r1;

        unsigned int croot1 = find_component_root(g, r1);
        unsigned int croot2 = find_component_root(g, r2);

        if (croot1 != croot2) {
            g->total_distance_div -= (g->component_sizes[croot1])*(g->component_sizes[croot1]-1);
            g->total_distance_div -= (g->component_sizes[croot2])*(g->component_sizes[croot2]-1);
            if (g->component_sizes[croot1] > g->component_sizes[croot2]) {
                g->component_tree[croot2] = croot1;
                g->component_sizes[croot1] += g->component_sizes[croot2];
                g->total_distance_div += (g->component_sizes[croot1])*(g->component_sizes[croot1]-1);
            } else {
                g->component_tree[croot1] = croot2;
                g->component_sizes[croot2] += g->component_sizes[croot1];
                g->total_distance_div += (g->component_sizes[croot2])*(g->component_sizes[croot2]-1);
            }
        }
        update_distance(g);
    }

    g->esize += 1;
}





void vertex_step(struct es_graph *g) {
    g->distribution[g->esize*2] = g->vsize;
    unsigned int r;
    if (rand_double() * (2*g->esize+1 + 2*g->d) >= 2*g->d) {
        r = g->distribution[rand_int(2*g->esize + 1)];
    } else {
        r = (unsigned int)rand_int(g->vsize+1);
    }
    struct edge *e = &(g->e[g->esize]);
    e->head = r;
    e->tail = g->vsize;
    g->distribution[2*g->esize+1] = r;
    struct node* newnode = init_node(&g->v[g->vsize]);
    newnode->parent = r;
    if (r != g->vsize) {
        newnode->nbs_size = 1;
        newnode->nbs_list[0] = r;
        unsigned int croot = find_component_root(g, r);
        g->component_tree[g->vsize] = croot;
        g->component_sizes[croot] += 1;
        g->v[r].nbs_size++;
        if (g->v[r].nbs_size > g->v[r].max_nbs_size) {
            resize_nbs(&g->v[r]);
        }
        g->v[r].nbs_list[g->v[r].nbs_size-1] = g->vsize;

        unsigned int prev = g->vsize;
        unsigned int tmp = r;
        unsigned int depth = 1;
        while (true) {
            if (depth+1 > g->v[tmp].max_tree_depth) {
                resize_tree(&g->v[tmp]);
            }

            if (depth + 1 > g->v[tmp].tree_depth) {
                g->v[tmp].tree_depth = depth + 1;
            }

            g->v[tmp].tree_struct[depth] += 1;
            g->v[tmp].num_children += 1;

            g->total_distance += depth;
            for (unsigned int i = 1; i < g->v[tmp].tree_depth; i++) {
                g->total_distance += (depth+i)*(g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i-1]);
            }

            if (g->v[tmp].parent == tmp) {
                break;
            }
            prev = tmp;
            tmp = g->v[tmp].parent;
            depth++;

        }

        bool* found = (bool*)malloc(sizeof(bool)*g->vsize);
        struct nd_tuple* q = malloc(sizeof(struct nd_tuple)*g->vsize);
        unsigned int qhead = 0;
        unsigned int qtail = 0;

        memset(found, false, sizeof(bool)*g->vsize);
        found[tmp] = true;
        q[qtail].node = tmp;
        q[qtail].dist = 0;
        qtail++;
        while (qhead != qtail) {
            struct nd_tuple cur = q[qhead];
            qhead++;


            if (cur.dist > 0) {
                for (unsigned int i = 0; i < g->v[cur.node].tree_depth; i++) {
                    g->total_distance += (depth + i + cur.dist) * g->v[cur.node].tree_struct[i];
                }
            }

            for (unsigned int j = 0; j < g->v[cur.node].nbs_size; j++) {
                unsigned int nb_id = g->v[cur.node].nbs_list[j];
                if (!found[nb_id] && g->v[nb_id].parent == nb_id) {
                    q[qtail].node = nb_id;
                    q[qtail].dist = cur.dist + 1;
                    qtail++;
                    found[nb_id] = true;
                }
            }
        }
        free(q);
        free(found);


        g->total_distance_div += 2*(g->component_sizes[croot]-1);


//        g->v[tmp].tree_struct[depth] += 1;
    } else {
        g->component_tree[g->vsize] = g->vsize;
        g->component_sizes[g->vsize] = 1;
    }
    g->vsize += 1;
    g->esize += 1;
}

void evolve(struct es_graph *g) {
    if (rand_double() < p()) {
        vertex_step(g);
    } else {
        edge_step(g);
    }
//    if (g->total_distance > 847208403) {
//        printf("here");
//
//        printf("%u\n%u\n", g->vsize, g->esize);
//        for (int i = 0; i < g->esize; i++) {
//            printf("(%u, %u) ", g->e[i].head, g->e[i].tail);
//        }
//        printf("\n");
//        for (int i = 0; i < g->vsize; i++) {
//            printf("\n%u (%u): ", i, g->v[i].nbs_size);
//            for (int j = 0; j < g->v[i].nbs_size; j++) {
//                printf("%u ", g->v[i].nbs_list[j]);
//            }
//        }
//        printf("\n");
//        for (int i = 0; i < g->vsize; i++) {
//            printf("\n%u (%s, %u, %u): ", i, (g->v[i].parent == i ? "T" : "F") , g->v[i].tree_depth, g->v[i].num_children);
//            for (int j = 0; j < g->v[i].tree_depth; j++) {
//                printf("%u ", g->v[i].tree_struct[j]);
//            }
//        }
//        printf("\n\n");
//        for (unsigned int i = 0; i < g->vsize; i++) {
//            if (g->component_tree[i] == i) {
//                printf("%u, ", g->component_sizes[i]);
//            }
//        }
//        printf("\n\nThe end\n\n\n");
//    }
}

unsigned long total_distance2(struct es_graph *g) {
    unsigned int* coeffs = malloc(sizeof(int)*g->vsize);
    memset(coeffs, 0, sizeof(int)*g->vsize);
    unsigned int* counters = malloc(sizeof(int)*g->vsize);
    bool* trimmed = calloc(g->vsize, sizeof(bool));
    bool* seen = malloc(sizeof(bool)*g->vsize);
    unsigned int* components = calloc(g->vsize, sizeof(unsigned int));
    unsigned int* path_lengths = calloc(g->vsize, sizeof(unsigned int));
    unsigned int diameter = 0;
    struct simple_doubly_linked_list *sll = malloc(sizeof(struct simple_doubly_linked_list));
    init_simple_linked_list(sll);

    for (unsigned int i = 0; i < g->vsize; i++) {
        counters[i] = g->v[i].nbs_size;
        coeffs[i] = 1;
    }

    // determine components
    unsigned int component_counter = 0;
    for (unsigned int i = 0; i < g->vsize; i++) {
        if (components[i] != 0) continue;

        component_counter++;
        simple_append(sll, i);
        unsigned int cur;
        while (!simple_is_empty(sll)) {
//            for (unsigned int k = 0; k < g->vsize; k++) {
//                printf("%3d: %d, %d, %d\n", k, trimmed[k], coeffs[k], counters[k]);
//            }
//            printf("\n\n");
            cur = simple_pop(sll);
            if (components[cur] != 0) {
                continue;
            }
            components[cur] = component_counter;
            for (unsigned int j = 0; j < g->v[cur].nbs_size; j++) {
                unsigned int nb_id = g->v[cur].nbs_list[j];
                if (components[nb_id] == 0) {
                    simple_append(sll, nb_id);
                }
            }
        }
    }

    unsigned int* component_sizes = calloc(component_counter+1, sizeof(unsigned int));
    for (unsigned int i = 0; i < g->vsize; i++) {
        component_sizes[components[i]]++;
    }
    for (unsigned int i = 0; i < g->vsize; i++) {
        components[i] = component_sizes[components[i]];
    }



    // Trim leaves
    for (unsigned int i = 0; i < g->vsize; i++) {
        if (trimmed[i] || counters[i] != 1) {
            continue;
        }

        simple_append(sll, i);
        unsigned int cur;
        while (!simple_is_empty(sll)) {
//            for (unsigned int k = 0; k < g->vsize; k++) {
//                printf("%3d: %d, %d, %d\n", k, trimmed[k], coeffs[k], counters[k]);
//            }
//            printf("\n\n");
            cur = simple_pop(sll);
            if (trimmed[cur]) {
                continue;
            }
            trimmed[cur] = true;
            for (unsigned int j = 0; j < g->v[cur].nbs_size; j++) {
                unsigned int nb_id = g->v[cur].nbs_list[j];
                if (trimmed[nb_id]) continue;
                counters[nb_id]--;
                coeffs[nb_id] += coeffs[cur];
                if (path_lengths[cur] + 1 > path_lengths[nb_id]) {
                    path_lengths[nb_id] = path_lengths[cur] + 1;
                }
                if (!trimmed[nb_id] && counters[nb_id] == 1) {
                    simple_append(sll, nb_id);
                }
            }
        }
    }

//    for (unsigned int k = 0; k < g->vsize; k++) {
//        printf("%3d: %d, %d, %d\n", k, trimmed[k], coeffs[k], counters[k]);
//    }
//    printf("\n\n");


    unsigned long total = 0;
    struct linked_list *ll = malloc(sizeof(struct linked_list));
    init_linked_list(ll);
    for (unsigned int i = 0; i < g->vsize; i++) {
        if (trimmed[i]) continue;
        memset(seen, false, sizeof(bool)*g->vsize);
        struct nd_tuple* n = malloc(sizeof(struct nd_tuple));
        n->node = i;
        n->dist = 0;
        append(ll, n);
        struct nd_tuple *cur;
        while ((cur = (struct nd_tuple*)pop(ll)) != NULL) {
            if (seen[cur->node] || trimmed[cur->node]) {
                free(cur);
                continue;
            }
//            printf("(%u, %u)\n", cur->node, cur->dist);
            seen[cur->node] = true;
            total += cur->dist * coeffs[cur->node]*coeffs[i];
            if (cur->dist + path_lengths[cur->node] + path_lengths[i] > diameter) {
                diameter = cur->dist + path_lengths[cur->node] + path_lengths[i];
            }
            for (unsigned int j = 0; j < g->v[cur->node].nbs_size; j++) {
                unsigned int nb_id = g->v[cur->node].nbs_list[j];
                if (!seen[nb_id] && !trimmed[nb_id]) {
                    struct nd_tuple *next = malloc(sizeof(struct nd_tuple));
                    next->node = nb_id;
                    next->dist = cur->dist + 1;
                    append(ll, next);
                }
            }
            free(cur);
        }

    }
    total >>= 1;

    for (unsigned int i = 0; i < g->vsize; i++) {
        if (!trimmed[i]) continue;
        total += coeffs[i]*(components[i]-coeffs[i]);
    }

    free(path_lengths);
    free(components);
    free(ll);
    free(sll);
    free(coeffs);
    free(counters);
    free(trimmed);
    free(seen);
    free(component_sizes);
    printf("%u\n", diameter);
    return total;
}

unsigned long total_distance(struct es_graph *g) {
    bool* seen = malloc(sizeof(bool)*g->vsize);
    unsigned long total = 0;
    struct linked_list *ll = malloc(sizeof(struct linked_list));
    init_linked_list(ll);
    for (unsigned int i = 0; i < g->vsize; i++) {
        memset(seen, false, sizeof(bool)*g->vsize);
        struct nd_tuple* n = malloc(sizeof(struct nd_tuple));
        n->node = i;
        n->dist = 0;
        append(ll, n);
        struct nd_tuple *cur;
        while ((cur = (struct nd_tuple*)pop(ll)) != NULL) {
            if (seen[cur->node]) {
                free(cur);
                continue;
            }
//            printf("(%u, %u)\n", cur->node, cur->dist);
            seen[cur->node] = true;
            total += cur->dist;
            for (unsigned int j = 0; j < g->v[cur->node].nbs_size; j++) {
                unsigned int nb_id = g->v[cur->node].nbs_list[j];
                if (!seen[nb_id]) {
                    struct nd_tuple *next = malloc(sizeof(struct nd_tuple));
                    next->node = nb_id;
                    next->dist = cur->dist + 1;
                    append(ll, next);
                }
            }
            free(cur);
        }

    }
    free(seen);
    free(ll);
    return total>>1;
}


double p() {
    return 0.90;
}

int main() {
//    printf("%f\n", 0.9);
//    printf("%d", false);
//    exit(0);
//    srandom(time(0));
//    srand48(time(0));
//    srand48(33112);
    srand48(96723);
    struct es_graph* g = generate_edge_step(0, 8000);
    printf("%u\n%u\n", g->vsize, g->esize);
//    for (int i = 0; i < g->esize; i++) {
//        printf("%u: (%u, %u),  ", i, g->e[i].head, g->e[i].tail);
//    }
//    printf("\n");
//    for (int i = 0; i < g->vsize; i++) {
//        printf("\n%u (%u): ", i, g->v[i].nbs_size);
//        for (int j = 0; j < g->v[i].nbs_size; j++) {
//            printf("%u ", g->v[i].nbs_list[j]);
//        }
//    }
//    printf("\n");
//    for (int i = 0; i < g->vsize; i++) {
//        printf("\n%u (%s, %u, %u): ", i, (g->v[i].parent == i ? "T" : "F") , g->v[i].tree_depth, g->v[i].num_children);
//        for (int j = 0; j < g->v[i].tree_depth; j++) {
//            printf("%u ", g->v[i].tree_struct[j]);
//        }
//    }
//    unsigned long div_fact = 0;
//    printf("\n\n");
//    for (unsigned int i = 0; i < g->vsize; i++) {
//        printf("%u: (%u, %u)\n",i, g->component_tree[i], g->component_sizes[i]);
//        if (g->component_tree[i] == i) {
//            div_fact += g->component_sizes[i]*(g->component_sizes[i]-1);
//        }
//    }
//    printf("\n\n");
////    printf("%lu\n", total_distance(g));
//    printf("%lu\n", total_distance2(g));
//    printf("%lu, %lu\n", div_fact, g->total_distance_div);

    fflush(stdout);
//    return 0;
    exit(0);
}

//int main() {
//    srandom(time(0));
//    long sum = 0;
//    int reps = 1000000000;
//    for (int i = 0; i < reps; i++) {
//        sum += (rand_double() < p() ? 1 : 0);
//    }
//    printf("%ld\n%.12f\n", sum, ((double)sum / (double)reps));
//    return 0;
//}
