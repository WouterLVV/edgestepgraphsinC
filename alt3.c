//
// Created by wouter on 29-10-19.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>


//// Default Values

#define NBS_INIT_SIZE 32
#define TREE_INIT_SIZE 32
#define MAX_DIAMETER 254

//// Flags

#define CSV_OUTPUT 1
#define DEBUG 0


//// Struct definitions

struct nd_tuple {
    unsigned int node;
    unsigned int dist;
};

struct node {
    unsigned int parent;
    unsigned int tree_depth;
    unsigned int max_tree_depth;
    unsigned int *tree_struct;
};

struct edge {
    unsigned int head;
    unsigned int tail;
};

struct component {
    unsigned int vsize;
    unsigned int esize;
    unsigned long *typical_distance;
    unsigned int diameter;
};

struct es_graph {
    double d;

    unsigned int vsize;
    unsigned int esize;

    struct node *v;
    struct edge *e;

    unsigned int *distribution;

    unsigned int num_components;
    unsigned int *component_list;
    struct component *components;
    unsigned int *component_tree;

    unsigned int dsize;
    unsigned int csize;
    unsigned int *central_map;
    unsigned int *central_list;
    unsigned char *distance;
};


//// Function definitions

struct es_graph *generate_edge_step(double d, unsigned int t, unsigned int e_central);

void evolve(struct es_graph *g);

double p(struct es_graph *g);

void vertex_step(struct es_graph *g);

void edge_step(struct es_graph *g);

void update_distance(struct es_graph *g);

void print_stats(struct es_graph *g);

void print_component_up(struct es_graph *g);

void edge_update_distances(struct es_graph *g, unsigned int n1, unsigned int n2, unsigned int croot);

bool edge_exists(unsigned int r1, unsigned int r2, struct es_graph *g);


//// Array resizers;

void resize_tree(struct node *n);


//// Standard functions

void resize_central(struct es_graph *g);

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
    n->max_tree_depth = TREE_INIT_SIZE;
    n->tree_depth = 1;
    n->tree_struct = calloc(n->max_tree_depth, sizeof(int));
    n->tree_struct[0] = 1;
    return n;
}

struct component *init_component(struct component *c) {
    c->vsize = 1;
    c->esize = 0;
    c->typical_distance = (unsigned long *) calloc(MAX_DIAMETER, sizeof(long));
    c->diameter = 0;
    return c;
}

void free_component(struct component *c) {
    free(c->typical_distance);
}

unsigned int merge_components(struct es_graph *g, unsigned int croot1, unsigned int croot2) {
    struct component *c1 = &g->components[croot1];
    struct component *c2 = &g->components[croot2];

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
    for (unsigned int i = 1; i <= c2->diameter; i++) {
        c1->typical_distance[i] += c2->typical_distance[i];
    }
    g->component_tree[croot2] = croot1;

    for (unsigned int i = 0; i < g->num_components; i++) {
        if (g->component_list[i] > croot2) {
            g->component_list[i - 1] = g->component_list[i];
        }
    }
    g->num_components -= 1;
    free_component(c2);
    return croot1;
}

unsigned int find_component_root(struct es_graph *g, unsigned int n) {
    while (g->component_tree[n] != n) {
        n = g->component_tree[n];
    }
    return n;
}

void create_new_component(struct es_graph *g, unsigned int root) {
    g->component_tree[root] = root;
    init_component(&g->components[root]);
    g->components[root].esize = 1;
    g->component_list[g->num_components] = root;
    g->num_components += 1;
}

unsigned char dist(struct es_graph *g, unsigned int n1, unsigned n2) {
    unsigned int c1 = g->central_map[n1];
    unsigned int c2 = g->central_map[n2];
    return g->distance[c1 * g->dsize + c2];
}

void setdist(struct es_graph *g, unsigned int n1, unsigned n2, unsigned char value) {
    unsigned int c1 = g->central_map[n1];
    unsigned int c2 = g->central_map[n2];
    g->distance[c1 * g->dsize + c2] = value;
    g->distance[c2 * g->dsize + c1] = value;
}

void resize_tree(struct node *n) {
    unsigned int old_size = n->max_tree_depth;
    n->max_tree_depth *= 2;
    n->tree_struct = realloc(n->tree_struct, sizeof(int) * n->max_tree_depth);
    memset(n->tree_struct + sizeof(int) * old_size, 0, sizeof(int) * old_size);
}






void init_es_graph(struct es_graph *g, double d, unsigned int t, unsigned int e_central) {
    g->d = d;

    g->v = (struct node *) calloc(t, sizeof(struct node));
    g->e = (struct edge *) malloc(sizeof(struct edge) * t);

    g->distribution = (unsigned int *) malloc(sizeof(int) * t * 2);

    g->num_components = 1;
    g->component_list = (unsigned int *) calloc(t, sizeof(int));
    g->components = (struct component *) malloc(sizeof(struct component) * t);
    g->component_tree = (unsigned int *) calloc(t, sizeof(int));

    g->central_map = calloc(t, sizeof(unsigned int));
    g->central_list = calloc(t, sizeof(unsigned int));
    g->distance = calloc((e_central * e_central), sizeof(char));

    g->vsize = 1;
    g->esize = 1;

    g->dsize = e_central;
    g->csize = 1;

    g->component_tree[0] = 0;
    g->component_list[0] = 0;
    init_component(&g->components[0]);

    g->e[0].head = 0;
    g->e[0].tail = 0;

    g->distribution[0] = 0;
    g->distribution[1] = 0;

    init_node(&g->v[0]);
}

struct es_graph *generate_edge_step(double d, unsigned int t, unsigned int e_central) {
    // Create and initialize all values of the struct
    struct es_graph *g = malloc(sizeof(struct es_graph));

    init_es_graph(g, d, t, e_central);

    // Evolve t times
    for (unsigned int i = 1; i < t; i++) {
        evolve(g);
        if (CSV_OUTPUT) {
            print_component_up(g);
        }
    }
    return g;
}

void print_component_up(struct es_graph *g) {
    unsigned long gtotal = 0;
    unsigned long gtotaldiv = 0;
    unsigned long ctotal = 0;
    unsigned long ctotaldiv = 0;
    unsigned int max_diam = 0;

    printf("%u", g->num_components);
    for (unsigned int i = 0; i < g->num_components; i++) {
        struct component *c = &g->components[g->component_list[i]];
        printf(";%u;%u;%u;%u", g->component_list[i], c->vsize, c->esize, c->diameter);
        ctotal = 0;
        for (unsigned int j = 1; j <= c->diameter; j++) {
            printf(";%lu", c->typical_distance[j]);
            ctotal += j*c->typical_distance[j];
        }
        ctotaldiv = (c->vsize * (c->vsize-1));
        printf(";%lu;%f", ctotal, (double) (2*ctotal)/ (double) ctotaldiv );

        gtotal += ctotal;
        gtotaldiv += ctotaldiv;

        if (c->diameter > max_diam) max_diam = c->diameter;
    }
    printf(";%u;%u;%lu;%f;%u", g->esize, g->vsize, gtotal, (double) (2 * gtotal) / (double) gtotaldiv, max_diam);
    printf("\n");
}



void evolve(struct es_graph *g) {
    if (rand_double() < p(g)) {
        vertex_step(g);
    } else {
        edge_step(g);
    }
}

bool edge_exists(unsigned int r1, unsigned int r2, struct es_graph *g) {
    if (r1 == r2) return true;
    if (g->v[r1].parent != r1 || g->v[r2].parent != r2) {
        return g->v[r1].parent == r2 || g->v[r2].parent == r1;
    }
    return dist(g, r1, r2) == 1;
}

struct nd_tuple vertex_tree_update(struct es_graph* g, struct component* c, unsigned int r) {
    struct nd_tuple t;
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
        c->typical_distance[depth] += 1;


        for (unsigned int i = 1; i < g->v[tmp].tree_depth; i++) {
            if (depth + i > c->diameter && g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1] > 0) {
                c->diameter = depth + i;
            }
            c->typical_distance[(depth + i)] += (g->v[tmp].tree_struct[i] - g->v[prev].tree_struct[i - 1]);
        }

        if (depth > c->diameter) c->diameter = depth;

        if (g->v[tmp].parent == tmp) break;

        prev = tmp;
        tmp = g->v[tmp].parent;
        depth++;
    }
    t.node = tmp;
    t.dist = depth;
    return t;
}

void vertex_central_update(struct es_graph *g, struct component *c, struct nd_tuple t) {
    unsigned int root_of_tree = t.node;
    unsigned int depth_of_new_node = t.dist;
    for (unsigned int i = 0; i < g->csize; i++) {
        unsigned int target = g->central_list[i];
        unsigned char distance = dist(g, root_of_tree, target);

        if (distance > 0) {

            if (depth_of_new_node + g->v[target].tree_depth - 1 + distance > c->diameter) {
                c->diameter = depth_of_new_node + g->v[target].tree_depth - 1 + distance;
            }

            for (unsigned int j = 0; j < g->v[target].tree_depth; j++) {
                c->typical_distance[depth_of_new_node + j + distance] += g->v[target].tree_struct[j];
            }

        }
    }
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
        unsigned int croot = find_component_root(g, r);
        g->component_tree[g->vsize] = croot;
        struct component *c = &g->components[croot];
        c->vsize += 1;
        c->esize += 1;

        struct nd_tuple t = vertex_tree_update(g, c, r);
        vertex_central_update(g, c, t);
    } else {
        create_new_component(g, g->vsize);
        g->central_map[g->vsize] = g->csize;
        g->central_list[g->csize] = g->vsize;
        g->csize++;
    }
    g->vsize += 1;
    g->esize += 1;
}


void add_vertex_to_central(struct es_graph* g, unsigned int n, unsigned int parent) {
    g->central_map[n] = g->csize;
    g->central_list[g->csize] = n;
    unsigned int n_index = g->csize;
    unsigned int p_index = g->central_map[parent];

    for (unsigned int i = 0; i < g->csize; i++) {
        g->distance[n_index*g->dsize + i] = g->distance[p_index*g->dsize + i] + 1;
        g->distance[i*g->dsize + n_index] = g->distance[i*g->dsize + p_index] + 1;
    }
    g->csize++;
    if (g->csize + 1 >= g->dsize) {
        resize_central(g);
    }
}

void resize_central(struct es_graph *g) {
    g->distance = realloc(g->distance, g->dsize * g->dsize * 4 * sizeof(char));

    for (signed int i = g->dsize-1; i >= 0; i--) {
        memmove(&g->distance[i*g->dsize*2], &g->distance[i*g->dsize], g->dsize * sizeof(char));
        memset(&g->distance[i*g->dsize*2+g->dsize], 0, g->dsize * sizeof(char));
    }
    memset(&g->distance[g->dsize*g->dsize*2], 0, g->dsize * g->dsize * 2 * sizeof(char));
    g->dsize *= 2;
//    printf("Resized!\n");
//    exit(1);
}


void un_treeify(struct es_graph *g, unsigned int node) {
    unsigned int stack[MAX_DIAMETER];
    unsigned int ptr = 0;

    stack[ptr] = node;
    ptr++;
    unsigned int tmp = node;
    while (g->v[tmp].parent != tmp) {
        unsigned int parent = g->v[tmp].parent;
        stack[ptr] = parent;
        ptr++;
        g->v[tmp].parent = tmp;
        tmp = parent;
    }
    ptr--;
    unsigned int prev = stack[ptr];

    while (ptr > 0) {
        ptr--;
        unsigned int next = stack[ptr];
        for (unsigned int i = 0; i < g->v[next].tree_depth; i++) {
            g->v[prev].tree_struct[i + 1] -= g->v[next].tree_struct[i];
        }
        while (g->v[prev].tree_struct[g->v[prev].tree_depth - 1] == 0) {
            g->v[prev].tree_depth--;
        }
        add_vertex_to_central(g, next, prev);
        prev = next;
    }
}

void edge_new_update(struct es_graph *g, struct edge *e) {
    unsigned int head = e->head;
    unsigned int tail = e->tail;

    un_treeify(g, head);
    un_treeify(g, tail);

    unsigned int croot1 = find_component_root(g, head);
    unsigned int croot2 = find_component_root(g, tail);

    if (croot1 != croot2) {
        croot1 = merge_components(g, croot1, croot2);
    }
    edge_update_distances(g, e->head, e->tail, croot1);
//        update_distance(g);
}

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

    bool exists = edge_exists(r1, r2, g);
    if (!exists) {
        edge_new_update(g, e);
    }

    g->esize += 1;
}


void edge_update_distances(struct es_graph *g, unsigned int n1, unsigned int n2, unsigned int croot) {
    unsigned int* marked1 = malloc(sizeof(int) * g->vsize);
    unsigned int* marked2 = malloc(sizeof(int) * g->vsize);

    struct component *c = &g->components[croot];


    marked1[0] = n2;
    marked2[0] = n1;
    unsigned int m1size = 1;
    unsigned int m2size = 1;


    unsigned int counttotal = 2;
    for (unsigned int i = 0; i < g->vsize; i++) {
        if (g->v[i].parent == i && i != n1 && i != n2) {
            counttotal++;
            if (dist(g, n1, i) > dist(g, n2, i) + 1) {
                marked1[m1size] = i;
                m1size++;
            }
            if (dist(g, n2, i) > dist(g, n1, i) + 1) {
                marked2[m2size] = i;
                m2size++;
            }
        }
    }

    for (unsigned int i = 0; i < m1size; i++) {
        for (unsigned int j = 0; j < m2size; j++) {
            unsigned int vi = marked1[i];
            unsigned int vj = marked2[j];
            unsigned char dij = dist(g, vi, vj);
            unsigned char newij = dist(g, vi, n2) + 1 + dist(g, n1, vj);

            if (dij == 0 || newij < dij) {
                setdist(g, vi, vj, newij);
                for (unsigned int k = 0; k < g->v[vi].tree_depth; k++) {
                    for (unsigned int l = 0; l < g->v[vj].tree_depth; l++) {
                        if (l + newij + k > c->diameter) c->diameter = l + newij + k;
                        if (dij > 0) c->typical_distance[l + dij + k] -= g->v[vi].tree_struct[k] * g->v[vj].tree_struct[l];
                        c->typical_distance[l + newij + k] += g->v[vi].tree_struct[k] * g->v[vj].tree_struct[l];
                    }
                }
            }
        }
    }

    while (c->typical_distance[c->diameter] == 0) {
        c->diameter--;
    }

    free(marked1);
    free(marked2);
}

double p(struct es_graph *g) {
    return 0.90;
}

//double p(struct es_graph *g) {
//    double t = (double)g->esize;
//    return 1./pow(t, 1.01);
//}
//
//double p(struct es_graph *g) {
//    double t = (double)g->esize;
//    return 1./log2(t);
//}

//double p(struct es_graph *g) {
//    double t = (double)g->esize;
//    return 1./log(t);
//}

int main() {
//    srand48(time(0));
//    srand48(33112);
//    long seed = time(0);
    long seed = 96723;
//    long seed = 96724;
    srand48(seed);
    struct es_graph *g = generate_edge_step(0, 10000, 100);
//    unsigned int ctr = 0;
//    for (unsigned int i = 0; i < g->vsize; i++) {
//        ctr += g->v[i].parent == i;
//    }
//    printf("%u\n", ctr);
//    printf("%u\n%u\n", g->vsize, g->esize);
//    for (int i = 0; i < g->esize; i++) {
//        printf("%u: (%u, %u),  ", i, g->e[i].head, g->e[i].tail);
//    }

    fflush(stdout);
    return 0;
}