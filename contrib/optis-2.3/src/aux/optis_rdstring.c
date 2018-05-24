# include <stdlib.h>

# include "../optis_global.h"

char* rd_string(FILE *file, int size)
{
    char c, cn, *str;
    int is_comment;

    str = (char*)malloc(size*sizeof(char));

    is_comment = 0;
    while(!feof(file))
    {
    	c = fgetc(file);

    	/* Se leu um espaco, tab, fim de linha descarta.*/
    	if( !( (c=='\n') || (c==9) || (c==32) ) )
    	{
    		if(is_comment == 0){
    			cn = fgetc(file);
    			/* Inicio do comentario*/
    			if( (c=='/') && ( cn == '*') ){
    				is_comment = 1;
    			}else{
    				fseek(file, -2, SEEK_CUR);
    				fscanf(file, "%s", str);
    				return str;
    			}

    		}else if(is_comment == 1 && c=='*'){
    			cn = fgetc(file);
    			if(cn == '/'){
    				is_comment = 0;
    			}else{
    				fseek(file, -1, SEEK_CUR);
    			}
    		}
    	}
    }

    return NULL;

}
