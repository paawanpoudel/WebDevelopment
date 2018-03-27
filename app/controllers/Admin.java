/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package controllers;

/**
 *
 * @author pp376
 */


import play.*;
import play.mvc.*;
import java.util.*;
import models.*;
import play.data.validation.Error;
import java.io.*;
import java.io.File;
import play.data.validation.*;
import notifiers.*;


@With(Secure.class)

public class Admin extends Controller {

    public static void index(){
         render();
    }

    public static void results(){
        render();
    }

     public static void query_genes(String gene_id){
        render();
    }

    public static void error_message(String message){
        render(message);
     }

    @Catch(IllegalStateException.class)
    public static void logIllegalState(Throwable throwable) {

        Logger.error("Illegal state %s…", throwable);
        error_message("missing values");
    }




    public static void createRiverPlots (String gene_List, String query_type, String plot_type, String rowv, String scale){


	String modified_gene_List=gene_List.replaceAll("[\\W*]","','").toUpperCase().replaceAll("\\''+",",").replaceAll("\\,+",",");

       	 // replace all the expecial charactcter
        Error gene_list_error = validation.required(gene_List).error;
        Error query_type_error = validation.required(query_type).error;
        Error plot_type_error = validation.required(plot_type).error;

         if(validation.hasErrors() ){
           params.flash(); // add http parameters to the flash scope
           validation.keep(); // keep the errors for the next request

           if( gene_list_error !=null) {
              error_message("Please insert either Ensembl or HGNC identifiers.");
           }

           if( query_type_error !=null){
              error_message("Please select either Gene or Trancript.");
           }
           if(plot_type_error!=null) {
              error_message("Please select either River plot or Heat map or Bar plot.");
           }

         error_message("Please insert either Ensembl or HGNC identifiers for your genes/transcripts and select some plotting option.");
          }

        create_R_plots plots=new create_R_plots();
        boolean result=plots.isRserveRunning();

        if(result){

                  System.out.println(modified_gene_List+"the modified gene_ist");
                  InputStream bais=plots.plot(modified_gene_List,plot_type, query_type, rowv, scale);

		  if(bais!=null){
                       response.current().contentType = "image/jpg";
                       renderBinary(bais);
                  }
                  if(bais==null){
                       error_message("Minimum 2 genes required to plot heat map of gene");
                  }
                  else{
                       error_message("Please insert valid Ensembl or HGNC identifiers.");
                  }
        }
        else {
             error_message("It is not possible to visualise data at the moment. Please try later.");
        }
    }
}


