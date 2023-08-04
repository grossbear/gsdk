#include <stdio.h>
#include <string.h>

bool conv_to_dec_data(char buffer[], float val)
{
    int int_val = *(int*)&val;
    char val_buf[64]={0};
    sprintf(val_buf,"%X",int_val);

    size_t len = strlen(val_buf);
    if(len > 8) return false;

    char prefix[] = {"0x"};
    char zeros_tab[] = {"00000000"};

    size_t pos = 8 - len;
    zeros_tab[pos] = '\0';

    strcpy(buffer, prefix);
    strcat(buffer, zeros_tab);
    strcat(buffer, val_buf);

    return true;
}


/*


    FILE * file = fopen("float_tab.h","wb");

    for(uint16t half = 0; half < 32768; half++)
    {
        float flt_val = m_half_to_float_ilm(half);
        char buffer[64]={0};
        conv_to_dec_data(buffer, flt_val);

        if( !(half%6) )
            fprintf(file,"\n");
        fprintf(file, "%s, ", buffer);
        //printf("flt_val = %s\n", buffer);
    }

    fclose(file);

*/
