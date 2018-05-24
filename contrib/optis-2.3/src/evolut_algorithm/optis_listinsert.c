# include <stdlib.h>

# include "../optis_global.h"

/* A custom doubly linked list implemenation */
/* Insert an element X into the list at location specified by NODE */
void insert (list *node, int x)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n List Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(EXIT_FAILURE);
    }
    temp = (list *)malloc(sizeof(list));
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}
