package models;

import java.io.*;
import java.awt.*;
import java.awt.Image;

// import for connecting to R
import org.rosuda.REngine.*;
import org.rosuda.REngine.Rserve.*;
import java.nio.charset.Charset;
import controllers.*;

public class create_R_plots extends Canvas {


     public static boolean isRserveRunning() {
          RConnection c=null;
          try {

                c = new RConnection();
              // 	c.close();
		return true;

		}
          catch (Exception e) {

               System.err.println("Error: First connect try failed with: "+e.getMessage());
               return false;
          }

          finally{
               if(c!=null){
                    try{
                         c.close();
                    }
                    catch(Exception e){
                         System.err.println("Error: Coulnd't clean up connection to R!");
                         e.printStackTrace();
                   }

          }
         // return false;


}
     }

// this method creates the plot in R

    public static InputStream plot (String gene_sel,String plot_type, String query_type, String rowv, String scale){

          String message="Cannot connect to R";
          InputStream ret=null;         //new ByteArrayInputStream(message.getBytes(Charset.forName("UTF-8")));

          //InputStream ret=null;
      //   ret = new FileInputStream("error.jpeg");
         RConnection c = null;
         try {

              c = new RConnection();
              String device="jpeg";
              String gene_sel1="c('"+gene_sel.toString()+"')";
              String gene_query="gene.sel<-"+gene_sel1;
              String arg=rowv+","+"'"+scale+"'";
              REXP xp = c.parseAndEval("try("+device+"('test.jpg',quality=100,res=300,width=10,height=10,units='in'))");
              c.assign(".tmp.", gene_query);
              REXP r = c.parseAndEval("try(eval(parse(text=.tmp.)),silent=TRUE)");
              if (r.inherits("try-error")) {

                   System.err.println("For the R Error: "+r.toString());
              }
              else {

                   //System.out.print("Success!!!!");


                   if(plot_type.matches("h") && query_type.matches("g")){

                        // heatmap cannot be creted for a single gene

                        String[] num_gene=gene_sel.split(",");
                        if(num_gene.length>1){
                           c.eval("create_heatmap_gene("+arg+" )");

                         }

                        else{
                             c.close();
                           //  String error_message="Error: Minimum 2 genes required to query heatmap gene";
                             return ret;

                        }


                   }

                   else if(plot_type.matches("h") && query_type.matches("t")){

                        c.eval("create_heatmap_transcript("+arg+")");
                   }


                  // creates the riverplot depending upon the user input
                   else if(plot_type.matches("r") && query_type.matches("g")){

                        c.eval("riverplot_gene()");

                   }

                   else if(plot_type.matches("r") && query_type.matches("t") ){

                        c.eval("riverplot_transcript()");

                   }

                    else if(plot_type.matches("b")){

                         // only take the first gene form the list if the user send mulitiple genes
                         String[] gene=gene_sel.split(",");
                         String gene_sel2 = gene[0].replaceAll("'","");
                         gene_query="gene.sel<-"+"\""+gene_sel2+"\"";


                         System.out.print(gene_query+"query for creating bar plot");
                         c.assign(".tmp.", gene_query);
                         r = c.parseAndEval("try(eval(parse(text=.tmp.)),silent=TRUE)");
                         c.eval("barplotIsoformExpression()");

                   }

              }
              //creates the heatmap, barplat, riverplot depending upon the user's input


              xp = c.parseAndEval("r=readBin('test.jpg','raw',1024*1024); unlink('test.jpg'); r");
              Image img = Toolkit.getDefaultToolkit().createImage(xp.asBytes());
              InputStream in=new ByteArrayInputStream(xp.asBytes());

	      // remove the variables to not pollute the workspace
              c.parseAndEval("rm(gene.sel)");
              c.parseAndEval("rm(index)");
              c.parseAndEval("rm(full.data.sel)");
              c.parseAndEval("rm(data2plot)");
              c.parseAndEval("rm(rowv)");
              c.parseAndEval("rm(cell.type.expression)");
              c.parseAndEval("rm(cum.sum)");
              c.parseAndEval("rm(max.fpkm.sum)");
              c.parseAndEval("rm(min.fpkm.sum)");
              c.parseAndEval("rm(white.line)");
              c.parseAndEval("rm(cell.type.coordinates)");
              c.parseAndEval("rm(mean.norm)");
              c.parseAndEval("rm(text)");
              c.parseAndEval("rm(mean.norm)");
              c.parseAndEval("rm(q)");
              c.parseAndEval("rm(trans.colours)");
              c.parseAndEval("rm(gene.expression)");
              c.parseAndEval("rm(sum.trans.values)");
              c.parseAndEval("rm(trans)");
              c.parseAndEval("rm(y.limit)");

              return in;


         }
         catch (RserveException rse) { // RserveException (transport layer - e.g. Rserve is not running)
     //
              System.out.println(rse);
              return ret;
         }
         catch (REXPMismatchException mme) { // REXP mismatch exception (we got something we didn't think we get)

              mme.printStackTrace();
              return ret;

        }
         catch(Exception e) { // something else

              System.out.println("Something went wrong, but it's not the Rserve: "+e.getMessage());
              e.printStackTrace();
              return ret;
        } finally{
              if(c != null){
                   try{
                        c.close();
                   } catch(Exception e){
                        System.err.println("Coulnd't clean up connection to R!");
                        e.printStackTrace();
                   }
              }
        }
    }



}
