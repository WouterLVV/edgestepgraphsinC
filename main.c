#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <string.h>


struct simple_linked_list_node {
    struct simple_linked_list_node* next;
    unsigned int data;
};

struct simple_linked_list {
    struct simple_linked_list_node* first;
    struct simple_linked_list_node* last;
};

void init_simple_linked_list (struct simple_linked_list* ll) {
    ll->first = NULL;
    ll->last = NULL;
};

void simple_append(struct simple_linked_list* sll, unsigned int data) {
    struct simple_linked_list_node *new_node = malloc(sizeof(struct simple_linked_list_node));
    new_node->data = data;
    new_node->next = NULL;
    if (sll->first == NULL) {
        sll->first = new_node;
        sll->last = new_node;
    } else {
        sll->last->next = new_node;
        sll->last = new_node;
    }
}

unsigned int simple_pop(struct simple_linked_list* ll) {
    if (ll->first == NULL) {
        return 0;
    }
    struct simple_linked_list_node *front = ll->first;
    unsigned int ret = front->data;
    ll->first = front->next;
    free(front);
    return ret;
}

bool simple_is_empty(struct simple_linked_list* ll) {
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
    unsigned int* nbs_list;
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
long rand_int(int n) {
    if ((n - 1) == RAND_MAX) {
        return lrand48();
    } else {
        // Supporting larger values for n would requires an even more
        // elaborate implementation that combines multiple calls to rand()
        assert (n <= RAND_MAX);

        // Chop off all of the values that would cause skew...
        int end = RAND_MAX / n; // truncate skew
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

unsigned long total_distance2(struct es_graph *g);

struct es_graph *generate_edge_step(double d, int t) {
    struct es_graph *g = malloc(sizeof(struct es_graph));
    g->vsize = 1;
    g->esize = 1;
    g->d = d;
    g->v = (struct node*)calloc(t, sizeof(struct node));
    g->e = (struct edge*)malloc(sizeof(struct edge)*t);
    g->distribution = (unsigned int*)malloc(sizeof(int)*t*2);

    g->e[0].head = 0;
    g->e[0].tail = 0;
    g->distribution[0] = 0;
    g->distribution[1] = 0;

    for (int i = 1; i < t; i++) {
        evolve(g);
    }
    build_nodes(g);
    collapse_neighbors(g);
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

double p() {
    return 0.90;
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
    g->vsize += 1;
    g->esize += 1;
}

void evolve(struct es_graph *g) {
    if (rand_double() < p()) {
        vertex_step(g);
    } else {
        edge_step(g);
    }
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
    struct simple_linked_list *sll = malloc(sizeof(struct simple_linked_list));
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



int main() {
//    printf("%f\n", 0.9);
//    printf("%d", false);
//    exit(0);
//    srandom(time(0));
    srand48(time(0));
//    srand48(333111);
    struct es_graph* g = generate_edge_step(0, 20000);
    printf("%u\n%u\n", g->vsize, g->esize);
//    for (int i = 0; i < g->esize; i++) {
//        printf("(%u, %u) ", g->e[i].head, g->e[i].tail);
//    }
//    printf("\n");
//    for (int i = 0; i < g->vsize; i++) {
//        printf("\n%u (%u): ", i, g->v[i].nbs_size);
//        for (int j = 0; j < g->v[i].nbs_size; j++) {
//            printf("%u ", g->v[i].nbs_list[j]);
//        }
//    }
//    printf("\n");
//    printf("%lu\n", total_distance(g));
    printf("%lu\n", total_distance2(g));
    fflush(stdout);
    return 0;
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
