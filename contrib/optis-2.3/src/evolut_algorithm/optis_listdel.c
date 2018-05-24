# include <stdlib.h>

# include "../optis_global.h"

/* Delete the node NODE from the list */
list* del (list *node)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n List Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(EXIT_FAILURE);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child!=NULL)
    {
        temp->child->parent = temp;
    }
    free (node);
    return (temp);
}
