package controllers;

import java.util.*;

import models.*;

public class Users extends Application {

 
    public static void show(Long id) {
        User user = User.findById(id);
        notFoundIfNull(user);
        render();
    }
}
