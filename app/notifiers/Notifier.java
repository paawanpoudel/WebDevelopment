/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor

/**
 *
 * @author pp376
 */
package notifiers;

import play.mvc.*;

import javax.mail.internet.*;

import models.*;

public class Notifier extends Mailer {

    public static boolean welcome(User user) throws Exception {
        setFrom(new InternetAddress("pp376@cam.ac.uk", "Administrator"));
        setReplyTo(new InternetAddress("pp376@cam.ac.uk", "Help"));
        setSubject("Welcome %s", user.email);
        addRecipient(user.email, new InternetAddress("pp376@cam.ac.uk", "New users notice"));
        return sendAndWait(user);
    }
    
}

