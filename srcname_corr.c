//correct source name in the header from nuppi data for dada format
//add B/J to the front, remove CAL_ at the front for calibrators

#include<string.h>

void srcname_corr(char *src_name, char *src_name_new)
{
  int len,i;
  char src_name_1713[10],src_name_0613[10],src_name_1012[10],src_name_1022[10],src_name_2317[10],src_name_1518[10],src_name_1643[10],src_name_0030[10], src_name_0751[10];

  strcpy(src_name_0751,"0751+18");
  strcpy(src_name_1713,"1713+07");
  strcpy(src_name_0613,"0613-01");
  strcpy(src_name_1012,"1012+53");
  strcpy(src_name_1022,"1022+10");
  strcpy(src_name_2317,"2317+14");
  strcpy(src_name_1518,"1518+49");
  strcpy(src_name_1643,"1643-12");
  strcpy(src_name_0030,"0030+04");

  len=strlen(src_name);

  //Convert from Bname to Jname
  if(strcmp(src_name,src_name_1713)==0) 
    {
      strcpy(src_name_new,"J1713+0747");
      return;
    }
  if(strcmp(src_name,src_name_0613)==0)
    {
      strcpy(src_name_new,"J0613-0200");
      return;
    }
  if(strcmp(src_name,src_name_1012)==0)
    {
      strcpy(src_name_new,"J1012+5307");
      return;
    }
  if(strcmp(src_name,src_name_1022)==0)
    {
      strcpy(src_name_new,"J1022+1001");
      return;
    }
  if(strcmp(src_name,src_name_2317)==0)
    {
      strcpy(src_name_new,"J2317+1439");
      return;
    }
  if(strcmp(src_name,src_name_1518)==0)
    {
      strcpy(src_name_new,"J1518+4904");
      return;
    }
  if(strcmp(src_name,src_name_1643)==0)
    {
      strcpy(src_name_new,"J1643-1224");
      return;
    }
  if(strcmp(src_name,src_name_0030)==0)
    {
      strcpy(src_name_new,"J0030+0451");
      return;
    }
  if(strcmp(src_name,src_name_0751)==0)
    {
      strcpy(src_name_new,"J0751+1807");
      return;
    }

  //Formularize source name
  //Add 'B' to Bname
  if(len==7 && strcmp(src_name,"3C309-1")!=0)
    {
      strcpy(src_name_new,"B");
      strcat(src_name_new,src_name);
    }
  //Add 'J' to Jname
  else if(len==9)
    {
      strcpy(src_name_new,"J");
      strcat(src_name_new,src_name);
    }
  //Deal with calibrators, with 'CAL_' before name
  else if(len==11)
    {
      src_name_new[0]='B';
      for(i=1;i<=7;i++) src_name_new[i]=src_name[i+3];
      src_name_new[8]='\0';
    }
  else if(len==13)
    {
      src_name_new[0]='J';
      for(i=1;i<=9;i++) src_name_new[i]=src_name[i+3];
      src_name_new[10]='\0';
    }
  else if(strcmp(src_name,"3C273")==0 || strcmp(src_name,"C273")==0)
    {
      strcpy(src_name_new,"3C273");
    }
  else if(strcmp(src_name,"3C309-1")==0)
    {
      strcpy(src_name_new,"3C309");
    }
  else
    {
      printf("Unrecognized source name format:%s!\n",src_name);
      exit(1);
    }
}
