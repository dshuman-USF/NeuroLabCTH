/*
   test vehicle for learning how to read excel .xls files
   and stuff results into some containers.
   NO support for .xlsx files.
   Thu Oct 30 14:40:10 EDT 2014
*/


#include "cth_cluster.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <freexl.h>


typedef vector <string> FileList;
typedef vector <FileList>::iterator VFileList;

int
main (int argc, char *argv[])
{
   int line_count = 0;
   unsigned int worksheet_index, page;
   unsigned short int idx;
   const void *handle = NULL;
   int ret;
   unsigned int max_worksheet;
   unsigned int rows;
   unsigned short columns;
   unsigned int row;
   unsigned short col;
   const char *worksheet_name;
   string    currline;
   bool      found = false;
   bool      error = false;
   stringstream sprd_data;
   bool told = false;
   string cell_text;
   size_t beg, end; 

   if (argc < 2)
   {
      fprintf(stderr, "usage: xlsrd file\n");
      return -1;
   }

   while (!error)
   {
         /* opening the .XLS file [Workbook] */
      ret = freexl_open(argv[1], &handle);
      printf("ret: %d\n  %p\n", ret,handle);
      if (ret != FREEXL_OK)
      {
         fprintf(stderr, "OPEN ERROR: %d\n", ret);
         return line_count;
      }

         /* querying BIFF Worksheet entries */
      ret = freexl_get_info (handle, FREEXL_BIFF_SHEET_COUNT, &max_worksheet);
      if (ret != FREEXL_OK)
      {
         fprintf (stderr, "GET-INFO [FREEXL_BIFF_SHEET_COUNT] Error: %d\n", ret);
         error = true;
         break;;
      }
      printf("Has %d pages\n",max_worksheet);
      
        // look for "all cells" or "Info File" page
      for (worksheet_index = 0; worksheet_index < max_worksheet && !error; worksheet_index++)
      {
         ret = freexl_get_worksheet_name (handle, worksheet_index, &worksheet_name);
printf("%s\n",worksheet_name);
         if (ret != FREEXL_OK)
         {
            fprintf (stderr, "GET-WORKSHEET-NAME Error: %d\n", ret);
            error = true;
            break;;
         }
         if (worksheet_name != NULL)
         { 
            if ( ((strstr(worksheet_name,FILE_XLS1) != NULL) ||
                 (strstr(worksheet_name,FILE_XLS2) != NULL)) ||
                (strstr(worksheet_name,INFO_XLS1) != NULL && 
                 strstr(worksheet_name,INFO_XLS2) == NULL))
            {
                 page = worksheet_index;
               found = true;
//               break;
            }
         }
         ret = freexl_select_active_worksheet(handle, worksheet_index);
         ret = freexl_worksheet_dimensions(handle, &rows, &columns);
         printf("Found rows: %d  cols: %d\n",rows,columns);
      }
      if (!found)
      {
         fprintf(stderr, "Unknown spreadsheet contents.  Only file list or info spreadsheets are handled by this program.\n");
         error = true;
         break;;
      }
         // select worksheet and get info
      ret = freexl_select_active_worksheet(handle, page);
      if (ret != FREEXL_OK)
      {
         fprintf (stderr, "SELECT-ACTIVE_WORKSHEET Error: %d\n", ret);
         error = true;
         break;;
      }
      ret = freexl_get_active_worksheet(handle, &idx);
      printf("Selected worksheet %d\n",idx);

      ret = freexl_worksheet_dimensions(handle, &rows, &columns);
      if (ret != FREEXL_OK)
      {
         fprintf (stderr, "WORKSHEET-DIMENSIONS Error: %d\n", ret);
         error = true;
         break;;
      }
      printf("Found rows: %d  cols: %d\n",rows,columns);

      freexl_get_worksheet_name (handle, page, &worksheet_name);
      printf("Worksheet name is: %s with  %d rows x %d columns\n", worksheet_name, rows, columns);

      for (row = 0; row < rows && ! error; row++)
      {
         FreeXL_CellValue cell;
         stringstream fields;
         vector<int> ftype(columns);
         fill(ftype.begin(),ftype.end(),0);
         currline.clear();
printf("\nrow %d:  \n",row);
         for (col = 0; col < columns && !error; col++) 
         {
            ret = freexl_get_cell_value(handle, row, col, &cell);
            if (ret != FREEXL_OK)
            {
               fprintf (stderr, "CELL-VALUE-ERROR (r=%u c=%u): %d\n", row, col, ret);
               error = true;
               break;;
            }
             // Read in the current line and keep track of the field types.
             // Later, we will apply some heuristics to see if the line
             // matches any signatures of lines we want to keep
            switch (cell.type)
            {
               case FREEXL_CELL_INT:
                   fields << cell.value.int_value << "\t";  
                   ftype[col] = FREEXL_CELL_INT;
printf("int\n");
                   break;
               case FREEXL_CELL_DOUBLE:
                   fields << setprecision(15) << cell.value.double_value << "\t";
                   ftype[col] = FREEXL_CELL_DOUBLE;
printf("double\n");
                   break;
               case FREEXL_CELL_TEXT:
               case FREEXL_CELL_SST_TEXT:
                   cell_text = cell.value.text_value;
cout << "[" << cell_text << "]" << endl;
                   beg = cell_text.find_first_not_of("' '\t");
                   end = cell_text.find_last_not_of("' '\t\n");
printf("%lu %lu\n",beg, end);
                   if (beg != string::npos)
                      cell_text=cell_text.substr(beg,end-beg+1);
                   else
                      cell_text.clear(); // empty or all ws
cout << cell_text << endl;
                   fields << cell_text  << "\t";
                   ftype[col] = FREEXL_CELL_TEXT;
printf("text\n");
                   break;
               case FREEXL_CELL_DATE:
               case FREEXL_CELL_DATETIME:
               case FREEXL_CELL_TIME:
printf("date\n");
                   fields << cell.value.text_value << "\t";  
                   ftype[col] = FREEXL_CELL_DATE;
                   break;
               case FREEXL_CELL_NULL:
               default:
                   fields << "\t";  
                   ftype[col] = FREEXL_CELL_NULL;
                   break;
            }
         }
         fields << endl;


         if (ftype[0] == FREEXL_CELL_TEXT && 
            (ftype[1] == FREEXL_CELL_INT || ftype[1] == FREEXL_CELL_DOUBLE) && 
             ftype[2] == FREEXL_CELL_DOUBLE && ftype[3] == FREEXL_CELL_DOUBLE)
         {
            sprd_data << fields.str();
            if (!told)
            {
               printf("An info spread\n");
               told = true;
            }
         }
         else if (ftype[0] == FREEXL_CELL_TEXT && fields.str()[0] == '/' &&
                  ftype[1] == FREEXL_CELL_TEXT && (ftype[2] == FREEXL_CELL_DOUBLE || ftype[2] == FREEXL_CELL_INT))
         {
            if (!told)
            {
               printf("A file spread\n");
               told = true;
            }
            sprd_data << fields.str();
            ++line_count;
         }
      }
      break;   // done
   }

cout << sprd_data.str();


   if (handle) /* closing the .XLS file [Workbook] */
   {
      ret = freexl_close(handle);
      if (ret != FREEXL_OK)
         fprintf (stderr, "CLOSE ERROR: %d\n", ret);
   }
   if (!error)
      return line_count;
   else
      return 0;
}

