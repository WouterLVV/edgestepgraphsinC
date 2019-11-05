//
// Created by wouter on 29-10-19.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <string.h>


//// Default Values

#define NBS_INIT_SIZE 32
#define TREE_INIT_SIZE 32
#define MAX_DIAMETER 1000

//// Flags

#define CSV_OUTPUT 1
#define DEBUG 0


//// Struct definitions

struct nd_tuple {
    unsigned int node;
    unsigned int dist;
};

struct node {
    unsigned int nbs_size;
    unsigned int max_nbs_size;
    unsigned int *nbs_list;
    unsigned int parent;
    unsigned int tree_depth;
    unsigned int max_tree_depth;
    unsigned int *tree_struct;
    unsigned int num_children;
};

struct edge {
    unsigned int head;
    unsigned int tail;
};

struct component {
    unsigned int vsize;
    unsigned int esize;
    unsigned long total_distance;
    unsigned long total_distance_div;
    unsigned long *typical_distance;
    unsigned int max_distance;
    unsigned int diameter;
};

struct es_graph {
    unsigned int vsize;
    unsigned int esize;
    double d;
    struct node *v;
    struct edge *e;
    unsigned int *distribution;
    unsigned int num_components;
    unsigned int *component_list;
    struct component *components;
//    unsigned int *component_sizes;
    unsigned int *component_tree;
    unsigned long total_distance;
    unsigned long total_distance_div;
    unsigned long *typical_distance;
    unsigned long max_distance;
    unsigned int diameter;
};


//// Standard functions

double rand_double() {
    return drand48();
}

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

struct node *init_node(struct node *n) {
    n->max_nbs_size = NBS_INIT_SIZE;
    n->max_tree_depth = TREE_INIT_SIZE;
    n->tree_depth = 1;
    n->nbs_size = 0;
    n->tree_struct = calloc(n->max_tree_depth, sizeof(int));
    n->tree_struct[0] = 1;
    n->nbs_list = malloc(sizeof(int) * n->max_nbs_size);
    n->num_children = 0;
    return n;
}

struct component *init_component(struct component *c) {
    c->vsize = 1;
    c->esize = 0;
    c->total_distance = 0;
    c->total_distance_div = 0;
    c->max_distance = MAX_DIAMETER;
    c->typical_distance = (unsigned long *) calloc(MAX_DIAMETER, sizeof(long));
    c->diameter = 0;
    return c;
}

void free_component(struct component *c) {
    free(c->typical_distance);
}

unsigned int merge_components(struct es_graph *g, unsigned int croot1, unsigned int croot2) {
    struct component* c1 = &g->components[croot1];
    struct component* c2 = &g->components[croot2];

    if (c2->diameter > c1->diameter) {
        struct component *tmpptr = c1;
        c1 = c2;
        c2 = tmpptr;
        unsigned int tmp = croot1;
        croot1 = croot2;
        croot2 = tmp;
    }

    c1->vsize += c2->vsize;
    c1->esize += c2->esize;
    c1->total_distance += c2->total_distance;
    g->total_distance_div -= (c1->total_distance_div + c2->total_distance_div);
    c1->total_distance_div = c1->vsize*(c1->vsize-1);
    g->total_distance_div += c1->total_distance_div;
    for (unsigned int i = 1; i <= c2->diameter; i++) {
        c1->typical_distance[i] += c2->typical_distance[i];
    }

    g->component_tree[croot2] = croot1;

    for (unsigned int i = 0; i < g->num_components; i++) {
        if (g->component_list[i] > croot2) {
            g->component_list[i-1] = g->component_list[i];
        }
    }
    g->num_components -= 1;
    free_component(c2);
    return croot1;
}



//// Function definitions

struct es_graph *generate_edge_step(double d, int t);

void evolve(struct es_graph *g);

double p();

void vertex_step(struct es_graph *g);

void edge_step(struct es_graph *g);

void update_distance(struct es_graph *g);

unsigned int find_component_root(struct es_graph *g, unsigned int n);

void print_stats(struct es_graph *g);

void mark_nodes(struct es_graph *g, unsigned int n1, unsigned int n2, unsigned int croot);


//// Array resizers;

void resize_tree(struct node *n);

void resize_typical(struct es_graph *g);

void resize_component_typical(struct component *c);

void resize_nbs(struct node *n);


struct es_graph *generate_edge_step(double d, int t) {
    // Create and initialize all values of the struct
    struct es_graph *g = malloc(sizeof(struct es_graph));
    g->vsize = 1;
    g->esize = 1;
    g->d = d;
    g->v = (struct node *) calloc(t, sizeof(struct node));
    g->e = (struct edge *) malloc(sizeof(struct edge) * t);
    g->distribution = (unsigned int *) malloc(sizeof(int) * t * 2);
//    g->component_sizes = (unsigned int *) malloc(sizeof(int) * t);
    g->num_components = 1;
    g->component_list = (unsigned int *) calloc(t, sizeof(int));
    g->components = (struct component *) malloc(sizeof(struct component)*t);
    g->component_tree = (unsigned int *) calloc(t, sizeof(int));
    g->total_distance = 0;
    g->total_distance_div = 0;
    g->diameter = 0;
    g->max_distance = MAX_DIAMETER;
    g->typical_distance = (unsigned long *) calloc(MAX_DIAMETER, sizeof(long));

//    g->component_sizes[0] = 1;
    g->component_tree[0] = 0;
    g->component_list[0] = 0;
    init_component(&g->components[0]);
    g->e[0].head = 0;
    g->e[0].tail = 0;
    g->distribution[0] = 0;
    g->distribution[1] = 0;
    init_node(&g->v[0]);

    // Evolve t times
    for (unsigned int i = 1; i < t; i++) {
        evolve(g);
        if (CSV_OUTPUT) {
            print_stats(g);
        }
    }
    return g;
}


void evolve(struct es_graph *g) {
    if (rand_double() < p()) {
        vertex_step(g);
    } else {
        edge_step(g);
    }
}

void print_stats(struct es_graph *g) {
    printf("%u;%lu;%f;%u", g->vsize, g->total_distance, (double) (2 * g->total_distance) / (double) g->total_distance_div,
           g->diameter);
    unsigned long sum = 0;
    for (unsigned int j = 1; j <= g->diameter; j++) {
        printf(";%lu", g->typical_distance[j]);
        sum += j*g->typical_distance[j];
    }
    printf(";%lu;;%u", sum, g->num_components);
    for (unsigned int i = 0; i < g->num_components; i++) {
        struct component *c = &g->components[g->component_list[i]];
        printf(";;%u;%lu;%f;%u", g->component_list[i], c->total_distance, (double) (2 * c->total_distance) / (double) c->total_distance_div,
               c->diameter);
    }
    printf("\n");
}


unsigned int find_component_root(struct es_graph *g, unsigned int n) {
    while (g->component_tree[n] != n) {
        n = g->component_tree[n];
    }
    return n;
}


//void update_distance(struct es_graph *g) {
//    bool *found = malloc(sizeof(bool) * g->vsize);
//    struct nd_tuple *q = malloc(sizeof(struct nd_tuple) * g->vsize);
//    unsigned int qhead = 0;
//    unsigned int qtail = 0;
//
//    // Reset typical distance and total
//    memset(g->typical_distance, 0, (g->diameter + 1) * sizeof(unsigned long));
//    unsigned long total = 0;
//
//    for (unsigned int i = 0; i < g->vsize; i++) {
//        if (g->v[i].parent != i) continue;
//        memset(found, false, sizeof(bool) * g->vsize);
//        found[i] = true;
//        qhead = 0;
//        qtail = 0;
//        q[qtail].node = i;
//        q[qtail].dist = 0;
//        qtail++;
//
//
//        while (qhead != qtail) {
//            struct nd_tuple cur = q[qhead];
//            qhead++;
//            total += cur.dist * (g->v[cur.node].num_children + 1) * (g->v[i].num_children + 1);
//
//            if (cur.dist + g->v[cur.node].tree_depth - 1 + g->v[i].tree_depth - 1 > g->diameter) {
//                g->diameter = cur.dist + g->v[cur.node].tree_depth - 1 + g->v[i].tree_depth - 1;
//                if (DEBUG) printf("370");
//                while (g->diameter + 1 > g->max_distance) {
//                    resize_typical(g);
//                }
//            }
//            if (cur.dist > 0) {
//                for (unsigned int j = 0; j < g->v[cur.node].tree_depth; j++) {
//                    for (unsigned int k = 0; k < g->v[i].tree_depth; k++) {
//                        g->typical_distance[k + cur.dist + j] += g->v[cur.node].tree_struct[j] * g->v[i].tree_struct[k];
//                    }
//                }
//            }
//
//            for (unsigned int j = 0; j < g->v[cur.node].nbs_size; j++) {
//                unsigned int nb_id = g->v[cur.node].nbs_list[j];
//                if (!found[nb_id] && g->v[nb_id].parent == nb_id) {
//                    q[qtail].node = nb_id;
//                    q[qtail].dist = cur.dist + 1;
//                    qtail++;
//                    found[nb_id] = true;
//                }
//            }
//        }
//
//    }
//    total >>= 1u;
//
//    for (unsigned int i = 0; i < g->vsize; i++) {
//        if (g->v[i].parent == i) continue;
//        total += (g->v[i].num_children + 1) *
//                 (g->component_sizes[find_component_root(g, i)] - (g->v[i].num_children + 1));
//        unsigned int parent = g->v[i].parent;
//        for (unsigned int j = 1; j < g->v[parent].tree_depth; j++) {
//            for (unsigned int k = 0; k < g->v[i].tree_depth; k++) {
//                if (k + j + 1 > g->diameter && g->v[parent].tree_struct[j] - g->v[i].tree_struct[j - 1] > 0) {
//                    g->diameter = k + j + 1;
//                    if (DEBUG) printf("405");
//                    while (g->diameter + 1 > g->max_distance) {
//                        resize_typical(g);
//                    }
//                }
//                g->typical_distance[(k + j + 1)] +=
//                        g->v[i].tree_struct[k] * (g->v[parent].tree_struct[j] - g->v[i].tree_struct[j - 1]);
//            }
//        }
//        for (unsigned int k = 0; k < g->v[i].tree_depth; k++) {
//            g->typical_distance[k + 1] += g->v[i].tree_struct[k];
//        }
//    }
//
//    for (unsigned int i = 0; i < g->vsize; i++) {
//        for (unsigned int j = 1; j < g->v[i].tree_depth; j++) {
//            g->typical_distance[j] += g->v[i].tree_struct[j];
//        }
//    }
//
//    for (unsigned int i = 1; i <= g->diameter; i++) {
//        g->typical_distance[i] >>= 1u;
//    }
//    free(q);
//    free(found);
//    g->total_distance = total;
//    while (g->typical_distance[g->diameter] == 0) {
//        g->diameter--;
//    }
//}


void edge_step(struct es_graph *g) {
    unsigned int r1, r2;
    if (rand_double() * (2 * g->esize + 2 * g->d) >= 2 * g->d) {
        r1 = g->distribution[rand_int(2 * g->esize)];
    } else {
        r1 = (unsigned int) rand_int(g->vsize);
    }
    if (rand_double() * (2 * g->esize + 2 * g->d) >= 2 * g->d) {
        r2 = g->distribution[rand_int(2 * g->esize)];
    } else {
        r2 = (unsigned int) rand_int(g->vsize);
    }
    struct edge *e = &g->e[g->esize];
    e->head = r1;
    e->tail = r2;
    g->distribution[2 * g->esize] = r1;
    g->distribution[2 * g->esize + 1] = r2;

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
        unsigned int q[g->diameter + 2];
        unsigned int qtail = 0;

        q[qtail] = r1;
        qtail++;
        unsigned int tmp = r1;
        while (g->v[tmp].parent != tmp) {
            unsigned int parent = g->v[tmp].parent;
            q[qtail] = parent;
            qtail++;
            g->v[tmp].parent = tmp;
            tmp = parent;
        }
        qtail--;
        unsigned int prev = q[qtail];

        while (qtail > 0) {
            qtail--;
            unsigned int next = q[qtail];
            for (unsigned int i = 0; i < g->v[next].tree_depth; i++) {
                g->v[prev].tree_struct[i + 1] -= g->v[next].tree_struct[i];
                g->v[prev].num_children -= g->v[next].tree_struct[i];
            }
            while (g->v[prev].tree_struct[g->v[prev].tree_depth - 1] == 0) {
                g->v[prev].tree_depth--;
            }
            prev = next;
        }
        q[qtail] = r2;
        qtail++;
        tmp = r2;
        while (g->v[tmp].parent != tmp) {
            unsigned int parent = g->v[tmp].parent;
            q[qtail] = parent;
            qtail++;
            g->v[tmp].parent = tmp;
            tmp = parent;
        }
        qtail--;
        prev = q[qtail];
        while (qtail > 0) {
            qtail--;
            unsigned int next = q[qtail];
            for (unsigned int i = 0; i < g->v[next].tree_depth; i++) {
                g->v[prev].tree_struct[i + 1] -= g->v[next].tree_struct[i];
                g->v[prev].num_children -= g->v[next].tree_struct[i];
            }
            while (g->v[prev].tree_struct[g->v[prev].tree_depth - 1] == 0) {
                g->v[prev].tree_depth--;
            }
            prev = next;
        }

        g->v[r1].nbs_size += 1;
        if (g->v[r1].nbs_size > g->v[r1].max_nbs_size) {
            resize_nbs(&g->v[r1]);
        }
        g->v[r1].nbs_list[g->v[r1].nbs_size - 1] = r2;

        g->v[r2].nbs_size += 1;
        if (g->v[r2].nbs_size > g->v[r2].max_nbs_size) {
            resize_nbs(&g->v[r2]);
        }
        g->v[r2].nbs_list[g->v[r2].nbs_size - 1] = r1;

        unsigned int croot1 = find_component_root(g, r1);
        unsigned int croot2 = find_component_root(g, r2);

        if (croot1 != croot2) {
            croot1 = merge_components(g, croot1, croot2);
        }
        mark_nodes(g, e->head, e->tail, croot1);
//        update_distance(g);

    }

    g->esize += 1;
}


void vertex_step(struct es_graph *g) {
    g->distribution[g->esize * 2] = g->vsize;
    unsigned int r;
    if (rand_double() * (2 * g->esize + 1 + 2 * g->d) >= 2 * g->d) {
        r = g->distribution[rand_int(2 * g->esize + 1)];
    } else {
        r = (unsigned int) rand_int(g->vsize + 1);
    }
    struct edge *e = &(g->e[g->esize]);
    e->head = r;
    e->tail = g->vsize;
    g->distribution[2 * g->esize + 1] = r;
    struct node *newnode = init_node(&g->v[g->vsize]);
    newnode->parent = r;
    if (r != g->vsize) {
        newnode->nbs_size = 1;
        newnode->nbs_list[0] = r;
        unsigned int croot = find_component_root(g, r);
        g->component_tree[g->vsize] = croot;
        struct component *c = &g->components[croot];
        c->vsize += 1;
        c->esize += 1;

        g->v[r].nbs_size++;
        if (g->v[r].nbs_size > g->v[r].max_nbs_size) {
            resize_nbs(&g->v[r]);
        }
        g->v[r].nbs_list[g->v[r].nbs_size - 1] = g->vsize;

        unsigned int prev = g->vsize;
        unsigned int tmp = r;
        unsigned int depth = 1;
        while (true) {
            if (depth + 1 > g->v[tmp].max_tree_depth) {
                resize_tree(&g->v[tmp]);
            }

            if (depth + 1 > g->v[tmp].tree_depth) {
                g->v[tmp].tree_depth = depth + 1;
            }

            g->v[tmp].tree_struct[depth] += 1;
            g->v[tmp].num_children += 1;

            g->total_distance += depth;
            g->typical_distance[depth] += 1;

            c->total_distance += depth;
            c->typical_distance[depth] += 1;


            for (unsigned int i = 1; i < g->v[tmp].tree_depth; i++) {
                g->total_distance += (depth + i) * (g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1]);
                c->total_distance += (depth + i) * (g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1]);
                if (depth + i > g->diameter && g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1] > 0) {
                    g->diameter = depth + i;
                    if (DEBUG) printf("598");
                    while (g->diameter + 1 > g->max_distance) {
                        resize_typical(g);
                    }
                }

                if (depth + i > c->diameter && g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1] > 0) {
                    c->diameter = depth + i;
                    if (DEBUG) printf("598");
                    while (c->diameter + 1 > c->max_distance) {
                        resize_component_typical(c);
                    }
                }

                g->typical_distance[(depth + i)] += (g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1]);
                c->typical_distance[(depth + i)] += (g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1]);

            }
            if (depth > g->diameter) {
                if (DEBUG) printf("607");
                g->diameter = depth;
                while (g->diameter + 1 > g->max_distance) {
                    resize_typical(g);
                }
            }

            if (depth > c->diameter) {
                if (DEBUG) printf("607");
                c->diameter = depth;
                while (c->diameter + 1 > c->max_distance) {
                    resize_component_typical(c);
                }
            }

            if (g->v[tmp].parent == tmp) {
                break;
            }
            prev = tmp;
            tmp = g->v[tmp].parent;
            depth++;

        }

        bool *found = (bool *) malloc(sizeof(bool) * g->vsize);
        struct nd_tuple *q = malloc(sizeof(struct nd_tuple) * g->vsize);
        unsigned int qhead = 0;
        unsigned int qtail = 0;

        memset(found, false, sizeof(bool) * g->vsize);
        found[tmp] = true;
        q[qtail].node = tmp;
        q[qtail].dist = 0;
        qtail++;
        while (qhead != qtail) {
            struct nd_tuple cur = q[qhead];
            qhead++;


            if (cur.dist > 0) {
                if (depth + g->v[cur.node].tree_depth - 1 + cur.dist > g->diameter) {
                    g->diameter = depth + g->v[cur.node].tree_depth - 1 + cur.dist;
                    if (DEBUG) printf("641");
                    while (g->diameter + 1 > g->max_distance) {
                        resize_typical(g);
                    }
                }

                if (depth + g->v[cur.node].tree_depth - 1 + cur.dist > c->diameter) {
                    c->diameter = depth + g->v[cur.node].tree_depth - 1 + cur.dist;
                    if (DEBUG) printf("641");
                    while (c->diameter + 1 > c->max_distance) {
                        resize_component_typical(c);
                    }
                }

                for (unsigned int i = 0; i < g->v[cur.node].tree_depth; i++) {
                    g->total_distance += (depth + i + cur.dist) * g->v[cur.node].tree_struct[i];
                    c->total_distance += (depth + i + cur.dist) * g->v[cur.node].tree_struct[i];

                    g->typical_distance[depth + i + cur.dist] += g->v[cur.node].tree_struct[i];
                    c->typical_distance[depth + i + cur.dist] += g->v[cur.node].tree_struct[i];
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


        g->total_distance_div += 2 * (c->vsize - 1);
        c->total_distance_div += 2 * (c->vsize - 1);


    } else {
        g->component_tree[g->vsize] = g->vsize;
        init_component(&g->components[g->vsize]);
        g->components[g->vsize].esize = 1;
        g->component_list[g->num_components] = g->vsize;
        g->num_components += 1;
    }
    g->vsize += 1;
    g->esize += 1;
}


void resize_tree(struct node *n) {
    unsigned int old_size = n->max_tree_depth;
    n->max_tree_depth *= 2;
    n->tree_struct = realloc(n->tree_struct, sizeof(int) * n->max_tree_depth);
    memset(n->tree_struct + sizeof(int) * old_size, 0, sizeof(int) * old_size);
}

void resize_typical(struct es_graph *g) {
    unsigned int old_size = g->max_distance;
    g->max_distance *= 2;
    g->typical_distance = realloc(g->typical_distance, sizeof(long) * g->max_distance);
    memset(g->typical_distance + sizeof(long) * old_size, 0, sizeof(long) * old_size);
}

void resize_component_typical(struct component *c) {
    unsigned int old_size = c->max_distance;
    c->max_distance *= 2;
    c->typical_distance = realloc(c->typical_distance, sizeof(long) * c->max_distance);
    memset(c->typical_distance + sizeof(long) * old_size, 0, sizeof(long) * old_size);
}


void resize_nbs(struct node *n) {
    n->max_nbs_size *= 2;
    n->nbs_list = realloc(n->nbs_list, sizeof(int) * n->max_nbs_size);
}

void mark_nodes(struct es_graph *g, unsigned int n1, unsigned int n2, unsigned int croot) {
    struct nd_tuple *q = malloc(sizeof(struct nd_tuple)*g->vsize);
    unsigned int* d1 = calloc(g->vsize, sizeof(unsigned int));
    unsigned int* d2 = calloc(g->vsize, sizeof(unsigned int));
    bool* found = (bool *) malloc(sizeof(bool) * g->vsize);
    unsigned int* marked1 = malloc(sizeof(int) * g->vsize);
    unsigned int* marked2 = malloc(sizeof(int) * g->vsize);
    bool* target = malloc(sizeof(bool)*g->vsize);
    memset(target, false, sizeof(bool)*g->vsize);
    struct nd_tuple cur;

    memset(found, false, sizeof(bool) * g->vsize);

    struct component *c = &g->components[croot];

    for (unsigned int i = 0; i < g-> vsize; i++) {
        d1[i] = -2;
        d2[i] = -2;
    }

//    memset(marked1, 255, sizeof(unsigned int) * g->vsize);
//    memset(marked2, 255, sizeof(unsigned int) * g->vsize);
    unsigned int qhead = 0;
    unsigned int qtail = 0;

    marked1[0] = n2;
    marked2[0] = n1;
    unsigned int m1size = 1;
    unsigned int m2size = 1;

    q[qtail].node = n1;
    q[qtail].dist = 0;
    qtail++;
    found[n1] = true;
    found[n2] = true;

    while (qhead != qtail) {
        cur = q[qhead];
        qhead++;
        d1[cur.node] = cur.dist;
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


    q[qtail].node = n2;
    q[qtail].dist = 0;
    qtail++;
    memset(found, false, sizeof(bool) * g->vsize);
    found[n2] = true;
    found[n1] = true;

    while (qhead != qtail) {
        cur = q[qhead];
        qhead++;
        d2[cur.node] = cur.dist;
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

    unsigned int counttotal = 2;
    for (unsigned int i = 0; i < g->vsize; i++) {
        if (g->v[i].parent == i && i != n1 && i != n2) {
            counttotal++;
            if (d1[i] > d2[i] + 1) {
                marked1[m1size] = i;
                m1size++;
            }
            if (d2[i] > d1[i] + 1) {
                marked2[m2size] = i;
                m2size++;
            }
        }
    }
//    printf("(%u, %u) of %u\n", m1size, m2size, counttotal);

    if (m2size < m1size) {
        unsigned int *tmpptr = marked2;
        unsigned int tmpsize = m2size;
        marked2 = marked1;
        m2size = m1size;
        marked1 = tmpptr;
        m1size = tmpsize;
    }


    for (unsigned int i = 0; i < m2size; i++) {
        target[marked2[i]] = true;
    }

    for (unsigned int i = 0; i < m1size; i++) {
        memset(found, false, sizeof(bool) * g->vsize);
        found[marked1[i]] = true;
        qhead = 0;
        qtail = 0;

        q[qtail].node = marked1[i];
        q[qtail].dist = 0;
        qtail++;

        while (qhead != qtail) {
            cur = q[qhead];
            qhead++;

            if (target[cur.node] && cur.dist > 0) {
//                g->total_distance -= cur.dist * (g->v[cur.node].num_children + 1) * (g->v[marked1[i]].num_children + 1);
                for (unsigned int j = 0; j < g->v[cur.node].tree_depth; j++) {
                    for (unsigned int k = 0; k < g->v[marked1[i]].tree_depth; k++) {
                        g->typical_distance[k + cur.dist + j] -= g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];
                        c->typical_distance[k + cur.dist + j] -= g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];

                        g->total_distance -= (k + cur.dist + j) * g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];
                        c->total_distance -= (k + cur.dist + j) * g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];
                    }
                }
            }

            for (unsigned int j = 0; j < g->v[cur.node].nbs_size; j++) {
                unsigned int nb_id = g->v[cur.node].nbs_list[j];
                if (!found[nb_id] && g->v[nb_id].parent == nb_id && !((cur.node == n1 && nb_id == n2) || (cur.node == n2 && nb_id == n1))) {
                    q[qtail].node = nb_id;
                    q[qtail].dist = cur.dist + 1;
                    qtail++;
                    found[nb_id] = true;
                }
            }
        }


        memset(found, false, sizeof(bool) * g->vsize);
        found[marked1[i]] = true;
        qhead = 0;
        qtail = 0;
        q[qtail].node = marked1[i];
        q[qtail].dist = 0;
        qtail++;


        while (qhead != qtail) {
            cur = q[qhead];
            qhead++;

            if (target[cur.node] && cur.dist > 0) {

                if (cur.dist + g->v[cur.node].tree_depth - 1 + g->v[marked1[i]].tree_depth - 1 > g->diameter) {
                    g->diameter = cur.dist + g->v[cur.node].tree_depth - 1 + g->v[marked1[i]].tree_depth - 1;
                    if (DEBUG) printf("722");
                    while (g->diameter + 1 > g->max_distance) {
                        resize_typical(g);
                    }
                }
                if (cur.dist + g->v[cur.node].tree_depth - 1 + g->v[marked1[i]].tree_depth - 1 > c->diameter) {
                    c->diameter = cur.dist + g->v[cur.node].tree_depth - 1 + g->v[marked1[i]].tree_depth - 1;
                    if (DEBUG) printf("722");
                    while (c->diameter + 1 > c->max_distance) {
                        resize_component_typical(c);
                    }
                }

                for (unsigned int j = 0; j < g->v[cur.node].tree_depth; j++) {
                    for (unsigned int k = 0; k < g->v[marked1[i]].tree_depth; k++) {
                        g->typical_distance[k + cur.dist + j] += g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];
                        c->typical_distance[k + cur.dist + j] += g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];

                        g->total_distance += (k + cur.dist + j)*g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];
                        c->total_distance += (k + cur.dist + j)*g->v[cur.node].tree_struct[j] * g->v[marked1[i]].tree_struct[k];
                    }
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
    }
    while (g->typical_distance[g->diameter] == 0) {
        g->diameter--;
    }
    while (c->typical_distance[c->diameter] == 0) {
        c->diameter--;
    }

    free(found);
    free(q);
    free(marked1);
    free(marked2);
    free(d1);
    free(d2);
    free(target);
}


double p() {
    return 0.90;
}

int main() {
//    srand48(time(0));
//    srand48(33112);
//    unsigned long seed = time(0);
    unsigned long seed = 96723;
    srand48(seed);
    struct es_graph *g = generate_edge_step(0, 40);
    printf("%u\n%u\n", g->vsize, g->esize);
    for (int i = 0; i < g->esize; i++) {
        printf("%u: (%u, %u),  ", i, g->e[i].head, g->e[i].tail);
    }
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
////    printf("%lu\n", total_distance2(g));
//    printf("%lu, %lu\n", div_fact, g->total_distance_div);

    fflush(stdout);
    return 0;
}