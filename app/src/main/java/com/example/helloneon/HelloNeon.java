package com.example.helloneon;

import java.io.InputStream;
import java.nio.ByteBuffer;

import android.content.res.AssetManager;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.widget.TextView;

public class HelloNeon extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_hello_neon);

        TextView textView = ((TextView) findViewById(R.id.text_view_hello_neon));
        //        textView.setText(stringFromJNI());

        try {
            byte[] argbBuffer = readBitmap();
            textView.setText(compareNeon(argbBuffer, 640, 480));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private byte[] readBitmap() throws Exception {
        AssetManager assetManager = getAssets();
        InputStream inputStream = assetManager.open("argb_640x480.jpeg");
        Bitmap bitmap = BitmapFactory.decodeStream(inputStream);
        int count = bitmap.getByteCount();
        ByteBuffer byteBuffer = ByteBuffer.allocate(count);
        bitmap.copyPixelsToBuffer(byteBuffer);

        return byteBuffer.array();
    }

    public native String stringFromJNI();

    public native String compareNeon(byte[] buffer, int width, int height);

    static {
        System.loadLibrary("hello-neon");
    }
}
