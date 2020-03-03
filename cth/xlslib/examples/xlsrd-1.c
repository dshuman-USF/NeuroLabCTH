/* 
/ xl2sql.c
/
/ FreeXL Sample code
/
/ Author: Sandro Furieri a.furieri@lqt.it
/
/ ------------------------------------------------------------------------------
/ 
/ Version: MPL 1.1/GPL 2.0/LGPL 2.1
/ 
/ The contents of this file are subject to the Mozilla Public License Version
/ 1.1 (the "License"); you may not use this file except in compliance with
/ the License. You may obtain a copy of the License at
/ http://www.mozilla.org/MPL/
/ 
/ Software distributed under the License is distributed on an "AS IS" basis,
/ WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
/ for the specific language governing rights and limitations under the
/ License.
/
/ The Original Code is the FreeXL library
/
/ The Initial Developer of the Original Code is Alessandro Furieri
/ 
/ Portions created by the Initial Developer are Copyright (C) 2011
/ the Initial Developer. All Rights Reserved.
/ 
/ Contributor(s):
/ Brad Hards
/ 
/ Alternatively, the contents of this file may be used under the terms of
/ either the GNU General Public License Version 2 or later (the "GPL"), or
/ the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
/ in which case the provisions of the GPL or the LGPL are applicable instead
/ of those above. If you wish to allow use of your version of this file only
/ under the terms of either the GPL or the LGPL, and not to allow others to
/ use your version of this file under the terms of the MPL, indicate your
/ decision by deleting the provisions above and replace them with the notice
/ and other provisions required by the GPL or the LGPL. If you do not delete
/ the provisions above, a recipient may use your version of this file under
/ the terms of any one of the MPL, the GPL or the LGPL.
/ 
*/


#include "../../cth_cluster.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "freexl.h"

int
main (int argc, char *argv[])
{
   unsigned int worksheet_index;
   const void *handle;
   int ret;
   unsigned int info;
   unsigned int max_worksheet;
   unsigned int rows;
   unsigned short columns;
   unsigned int row;
   unsigned short col;
   const char *utf8_string;
   const char *our_page="all cells";

   if (argc != 2)
   {
      fprintf(stderr, "usage: xlsrd path.xls\n");
      return -1;
   }

     /* opening the .XLS file [Workbook] */
   ret = freexl_open (argv[1], &handle);
   if (ret != FREEXL_OK)
   {
      fprintf(stderr, "OPEN ERROR: %d\n", ret);
      return -1;
   }

      /* checking for Password (obfuscated/encrypted) */
   ret = freexl_get_info (handle, FREEXL_BIFF_PASSWORD, &info);
   if (ret != FREEXL_OK)
   {
       fprintf (stderr, "GET-INFO [FREEXL_BIFF_PASSWORD] Error: %d\n", ret);
       goto stop;
   }
   switch (info)
   {
      case FREEXL_BIFF_PLAIN:
         break;
      case FREEXL_BIFF_OBFUSCATED:
      default:
         fprintf (stderr, "Password protected: (not accessible)\n");
         goto stop;
         break;
   }

      /* querying BIFF Worksheet entries */
   ret = freexl_get_info (handle, FREEXL_BIFF_SHEET_COUNT, &max_worksheet);
   if (ret != FREEXL_OK)
   {
      fprintf (stderr, "GET-INFO [FREEXL_BIFF_SHEET_COUNT] Error: %d\n", ret);
      goto stop;
   }

   for (worksheet_index = 0; worksheet_index < max_worksheet; worksheet_index++)
   {
      ret = freexl_get_worksheet_name (handle, worksheet_index, &utf8_string);
      if (ret != FREEXL_OK)
      {
         fprintf (stderr, "GET-WORKSHEET-NAME Error: %d\n", ret);
         goto stop;
      }
      if (utf8_string != NULL)  // we only want the "all cells" worksheet if the data_files spread
      { 
//       if (strstr(utf8_string,our_page) != NULL)
            printf ("%3u] WORKSHEET NAME: %s\n", worksheet_index, utf8_string);
//       else
//          continue;
      }
      else
         continue;

     /* selecting the active Worksheet */
      ret = freexl_select_active_worksheet (handle, worksheet_index);
      if (ret != FREEXL_OK)
      {
         fprintf (stderr, "SELECT-ACTIVE_WORKSHEET Error: %d\n", ret);
         goto stop;
      }
        /* dimensions */
      ret = freexl_worksheet_dimensions (handle, &rows, &columns);
      if (ret != FREEXL_OK)
      {
         fprintf (stderr, "WORKSHEET-DIMENSIONS Error: %d\n", ret);
         goto stop;
      }
      printf("Worksheet is %d rows x %d columns\n", rows, columns);

      for (row = 0; row < rows; row++)
      {
         FreeXL_CellValue cell;
//         for (col = 0; col < columns; col++)
         for (col = 0; col < 5; col++) // only want 1st 5
         {
            ret = freexl_get_cell_value (handle, row, col, &cell);
            if (ret != FREEXL_OK)
            {
               fprintf (stderr, "CELL-VALUE-ERROR (r=%u c=%u): %d\n", row, col, ret);
               goto stop;
            }
            switch (cell.type)
            {
               case FREEXL_CELL_INT:
                   printf ("%d", cell.value.int_value);
                   break;
               case FREEXL_CELL_DOUBLE:
                   printf ("%1.12f", cell.value.double_value);
                   break;
               case FREEXL_CELL_TEXT:
               case FREEXL_CELL_SST_TEXT:
                   printf("%s",cell.value.text_value);
                   break;
               case FREEXL_CELL_DATE:
               case FREEXL_CELL_DATETIME:
               case FREEXL_CELL_TIME:
                   printf ("'%s'", cell.value.text_value);
                   break;
               case FREEXL_CELL_NULL:
               default:
//                   printf (",");
                   break;

            }
            printf ("\t");
         }
         printf ("\n");
      }
   }


  stop:
/* closing the .XLS file [Workbook] */
    ret = freexl_close (handle);
    if (ret != FREEXL_OK)
      {
     fprintf (stderr, "CLOSE ERROR: %d\n", ret);
     return -1;
      }

    return 0;
}












