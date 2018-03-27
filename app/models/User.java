package models;

import javax.persistence.*;
import java.util.*;

import play.*;
import play.db.jpa.*;
import play.libs.*;
import play.data.validation.*;

@Entity
public class User extends Model {

    @Email
    @Required
    public String email;
    
    @Required
    public String passwordHash;
    
//    @Required
 //   public String name;
    
    public boolean isAdmin;
    
    public String needConfirmation;
    
    // ~~~~~~~~~~~~ 
    
    public User(String email, String password) {
        this.email = email;
        this.passwordHash = Codec.hexMD5(password);
    //    this.name = name;
        this.needConfirmation = Codec.UUID();
        create();
    }
    
    // ~~~~~~~~~~~~ 
    
    public boolean checkPassword(String password) {
        return passwordHash.equals(Codec.hexMD5(password));
    }

    public boolean isAdmin() {
        return email.equals(Play.configuration.getProperty("forum.adminEmail", ""));
    }
    
    // ~~~~~~~~~~~~ 
   
    public static User connect(String email, String password) {
        String password_hash=Codec.hexMD5(password);
       // System.out.println("the password is "+password+" the hashed password is"+password_hash);
        return find("byEmailAndPasswordHash", email, password_hash).first();
    }
        
    
    public static User findByEmail(String email) {
        return find("email", email).first();
    }

    public static User findByRegistrationUUID(String uuid) {
        return find("needConfirmation", uuid).first();
    }

    public static List<User> findAll(int page, int pageSize) {
        return User.all().fetch(page, pageSize);
    }

    public static boolean isEmailAvailable(String email) {
        return findByEmail(email) == null;
    }
    
}
