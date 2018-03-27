package controllers;


import play.*;
import play.mvc.*;
import java.util.*;
import models.*;
import play.data.validation.Error;
import java.io.*;
import java.io.File;
import play.data.validation.*;
import notifiers.*;

        
public class Application extends Controller {

   
      @Before
    static void globals() {
        renderArgs.put("connected", connectedUser());
    }

    public static void signup() {
        render();
    }

// added pages       
    public static void contact(){
         render();
    }

    public static void cite_us(){
        render();
    }
    public static void funding(){
        render();
     }
    
/// end
       

    public static void register(@Required @Email String email, @Required @MinSize(8) String password, @Equals("password") String password2) {
        if (validation.hasErrors()) {
            validation.keep();
            params.flash();
            flash.error("Please correct these errors !");
            signup();
        }
        User user = new User(email, password);
        try {
            if (Notifier.welcome(user)) {
                flash.success("Your account is created. Please check your emails ...");
                login();
            }
        } catch (Exception e) {
            Logger.error(e, "Mail error");
        }
        flash.error("Oops ... (the email cannot be sent)");
        login();
    }

    public static void confirmRegistration(String uuid) {
        User user = User.findByRegistrationUUID(uuid);
        notFoundIfNull(user);
        user.needConfirmation = null;
        user.save();
        connect(user);
        flash.success("Welcome %s !", user.email);
        Admin.index();
    }
    public static void login() {
        render();
    }

    public static void index(){
	Application.login();
    }    

    public static void logout() {
        flash.success("You've been logged out");
        session.clear();
        Application.login();
    }

    public static void resendConfirmation(String uuid) {
        User user = User.findByRegistrationUUID(uuid);
        notFoundIfNull(user);
        try {
            if (Notifier.welcome(user)) {
                flash.success("Please check your emails ...");
                flash.put("email", user.email);
                login();
            }
        } catch (Exception e) {
            Logger.error(e, "Mail error");
        }
        flash.error("Oops (the email cannot be sent)...");
        flash.put("email", user.email);
        login();
    }
    
    // ~~~~~~~~~~~~ Some utils
    
   
    static void connect(User user) {
        session.put("logged", user.id);
    }

    static User connectedUser() {
        String userId = session.get("logged");
        
        try{
             return userId == null ? null : (User) User.findById(Long.parseLong(userId));
        }
        catch (Exception e){
             return null;
        }
       
        
    }
    
    
}
