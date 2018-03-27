package controllers;

import java.lang.reflect.InvocationTargetException;
import java.util.List;
import java.util.Date;
import play.Play;
import play.mvc.*;
import play.data.validation.*;
import play.libs.*;
import play.utils.*;

import models.*;

public class Security extends Secure.Security {

    static boolean authentify(String username, String password) {
        return User.connect(username, password) != null;
    }
  
    static boolean check(String profile) {
   	 if("admin".equals(profile)) {
         	return User.find("byEmail", connected()).<User>first().isAdmin;

    	}
    return false;

   }

   static void onAuthenticated() {

        Admin.index();
   }
      
   static void onDisconnected() {
        Application.login();
   }
    
}
   

