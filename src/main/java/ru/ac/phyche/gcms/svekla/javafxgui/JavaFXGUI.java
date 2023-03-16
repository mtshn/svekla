package ru.ac.phyche.gcms.svekla.javafxgui;

import java.io.File;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.application.Application;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.concurrent.Worker;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.web.WebView;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.util.Duration;
import javafx.fxml.FXMLLoader;
import netscape.javascript.JSObject;

public class JavaFXGUI extends Application {

	@Override
	public void start(Stage primaryStage) throws Exception {
		primaryStage.setTitle("SVEKLA: A GUI for GC-MS retention index and mass spectra prediction");
		FXMLLoader loader = new FXMLLoader();
		loader.setLocation(new File("./svekla.fxml").toURI().toURL());
		VBox vBox = loader.<VBox>load();
		Scene scene = new Scene(vBox, 800, 1000);
		primaryStage.setScene(scene);
		primaryStage.show();
		primaryStage.setOnCloseRequest(event -> {
			System.exit(0);
		});

		TextField accessResults = (TextField) vBox.lookup("#accessresults");
		TextField nonpolarResult = (TextField) vBox.lookup("#rinonpolar");
		TextField sveklaResults = (TextField) vBox.lookup("#sveklaresults");
		TextField polarResult = (TextField) vBox.lookup("#ripolar");
		TextField nonpolarDev = (TextField) vBox.lookup("#devnonpolar");
		TextField polarDev = (TextField) vBox.lookup("#devpolar");
		TextField nonpolarTarget = (TextField) vBox.lookup("#targetnonpolar");
		TextField polarTarget = (TextField) vBox.lookup("#targetpolar");
		TextField spectralSimilarity = (TextField) vBox.lookup("#spectralsimilarity");
		TextField fileWithSpectrum = (TextField) vBox.lookup("#filespectrum");
		TextField molecularFormula = (TextField) vBox.lookup("#molecularformula");
		TextField targetMW = (TextField) vBox.lookup("#targetmw");

		TextArea targetMS = (TextArea) vBox.lookup("#textareaspectrum");
		TextArea textareaMSMS = (TextArea) vBox.lookup("#msmstextarea");
		CheckBox checkbox1 = (CheckBox) vBox.lookup("#checkboxsquareroot");
		CheckBox checkbox2 = (CheckBox) vBox.lookup("#checkboxallpeaks");
		CheckBox checkboxMSMS1 = (CheckBox) vBox.lookup("#checkboxmsmspos");
		CheckBox checkboxMSMS2 = (CheckBox) vBox.lookup("#checkboxmsmsneg");
		CheckBox checkboxMSMS3 = (CheckBox) vBox.lookup("#checkboxmsmsnewwin");
		CheckBox checkboxMSMS4 = (CheckBox) vBox.lookup("#checkboxmsmslargest");

		WebView webView = (WebView) vBox.lookup("#webview");
		Canvas canvas1 = (Canvas) vBox.lookup("#canvas1");
		Canvas canvas2 = (Canvas) vBox.lookup("#canvas2");
		Button button = (Button) vBox.lookup("#button");
		Button button2 = (Button) vBox.lookup("#button2");
		Button buttonLoad = (Button) vBox.lookup("#loadbutton");
		Button buttonMSMS = (Button) vBox.lookup("#buttonmsms");
		Button buttonUseSpectrum = (Button) vBox.lookup("#usepastebutton");
		Button buttonCloseSpectra = (Button) vBox.lookup("#closespectra");

		GraphicsContext g = canvas1.getGraphicsContext2D();
		g.setFill(Color.valueOf("#FFFFFF"));
		g.fillRect(1, 1, canvas1.getWidth() - 1, canvas1.getWidth() - 1);
		g = canvas2.getGraphicsContext2D();
		g.setFill(Color.valueOf("#FFFFFF"));
		g.fillRect(1, 1, canvas2.getWidth() - 1, canvas2.getWidth() - 1);

		String path = new java.io.File(".").getCanonicalPath();
		webView.getEngine().load("file://" + path + "/moledit.html");
		final JSJavaCall jsJavaCall = new JSJavaCall();
		webView.getEngine().getLoadWorker().stateProperty().addListener(new ChangeListener<Worker.State>() {
			@Override
			public void changed(@SuppressWarnings("rawtypes") ObservableValue observable, Worker.State oldValue,
					Worker.State newValue) {
				if (newValue != Worker.State.SUCCEEDED) {
					return;
				}
				JSObject window = (JSObject) webView.getEngine().executeScript("window");
				window.setMember("o", jsJavaCall);
			}
		});

		MSpectra.SpectrumMS[] target = new MSpectra.SpectrumMS[1];

		EventHandler<ActionEvent> predictAndRefreshSpectra = (new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				MSpectra.SpectrumMS sp = MSpectra.sp(jsJavaCall.smiles);
				int min = sp.minMZ;
				int max = sp.maxMZ;
				MSpectra.SpectrumMS comp = new MSpectra.SpectrumMS();
				if (target[0] != null) {
					min = Math.max(sp.minMZ, target[0].minMZ);
					max = Math.max(sp.maxMZ, target[0].maxMZ);
					MSpectra.displaySpectrumOnCanvas(canvas1, target[0], min, max, sp,
							checkbox1.selectedProperty().get(), checkbox2.selectedProperty().get(), "Unknown analyte");
					spectralSimilarity.setText(MSpectra.spectralSimilarity(sp, target[0]));
					comp = target[0];
				}
				String caption = (jsJavaCall.smiles.length() < 30) && (jsJavaCall.smiles.length() > 3)
						? "Predicted for SMILES:" + jsJavaCall.smiles
						: "Predicted";
				MSpectra.displaySpectrumOnCanvas(canvas2, sp, min, max, comp, checkbox1.selectedProperty().get(),
						checkbox2.selectedProperty().get(), caption);
			}
		});

		button.setOnAction(predictAndRefreshSpectra);
		checkbox1.setOnAction(predictAndRefreshSpectra);
		checkbox2.setOnAction(predictAndRefreshSpectra);

		button2.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				FileChooser fileChooser = new FileChooser();
				File f = fileChooser.showOpenDialog(primaryStage);
				fileWithSpectrum.setText(f.getAbsolutePath());
			}
		});
		buttonLoad.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				target[0] = MSpectra.load(fileWithSpectrum.getText().trim());
				if (target[0] != null) {
					int mw = 0;
					try {
						mw = (int) Float.parseFloat(targetMW.getText());
					} catch (Throwable e) {
					}
					if (mw > 0) {
						target[0].maxMZ = Math.min(target[0].maxMZ, mw + 10);
					}
					MSpectra.displaySpectrumOnCanvas(canvas1, target[0], target[0].minMZ, target[0].maxMZ,
							new MSpectra.SpectrumMS(), checkbox1.selectedProperty().get(),
							checkbox2.selectedProperty().get(), "Unknown analyte");
				}
			}
		});
		buttonUseSpectrum.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String text = targetMS.getText().replaceAll("\n", System.getProperty("line.separator"));
				String[] lines = text.split("\\r?\\n|\\)\\s+\\(|\\)\\(");
				target[0] = MSpectra.load(lines);
				if (target[0] != null) {
					int mw = 0;
					try {
						mw = (int) Float.parseFloat(targetMW.getText());
					} catch (Throwable e) {
					}
					if (mw > 0) {
						target[0].maxMZ = Math.min(target[0].maxMZ, mw + 10);
					}
					MSpectra.displaySpectrumOnCanvas(canvas1, target[0], target[0].minMZ, target[0].maxMZ,
							new MSpectra.SpectrumMS(), checkbox1.selectedProperty().get(),
							checkbox2.selectedProperty().get(), "Unknown analyte");
				}
			}
		});

		buttonMSMS.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String msms = "";
				if (checkboxMSMS1.selectedProperty().get()) {
					msms = msms + MSMSCFM.sp(jsJavaCall.smiles, true, checkboxMSMS4.selectedProperty().get());
				}
				if (checkboxMSMS2.selectedProperty().get()) {
					msms = msms + MSMSCFM.sp(jsJavaCall.smiles, false, checkboxMSMS4.selectedProperty().get());
				}
				textareaMSMS.textProperty().set(msms);
				if (checkboxMSMS3.selectedProperty().get()) {
					Stage stage = new Stage();
					stage.setTitle("Predicted MS/MS spectrum");
					TextArea ta = new TextArea();
					BorderPane brd = new BorderPane();
					stage.setScene(new Scene(brd, 600, 800));
					stage.sizeToScene();
					brd.setCenter(ta);
					ta.setText(msms);
					stage.show();
				}
			}
		});

		buttonCloseSpectra.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				GraphicsContext g = canvas1.getGraphicsContext2D();
				g.setFill(Color.valueOf("#FFFFFF"));
				g.fillRect(1, 1, canvas1.getWidth() - 1, canvas1.getWidth() - 1);
				target[0] = null;
			}
		});

		primaryStage.show();

		RIPrediction rp = new RIPrediction();

		nonpolarTarget.textProperty().addListener((observable, oldValue, newValue) -> {
			jsJavaCall.oldSmiles = "xxx";
		});

		polarTarget.textProperty().addListener((observable, oldValue, newValue) -> {
			jsJavaCall.oldSmiles = "xxx";
		});

		Timeline riTimer = new Timeline(new KeyFrame(Duration.seconds(0.5), new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				if (!jsJavaCall.smiles.equals(jsJavaCall.oldSmiles)) {
					rp.predictAndWrite(jsJavaCall.smiles, accessResults, nonpolarResult, sveklaResults, polarResult,
							nonpolarDev, polarDev, nonpolarTarget, polarTarget, molecularFormula);
					jsJavaCall.oldSmiles = jsJavaCall.smiles;
					GraphicsContext g = canvas2.getGraphicsContext2D();
					g.setFill(Color.valueOf("#FFFFFF"));
					g.fillRect(1, 1, canvas2.getWidth() - 1, canvas2.getWidth() - 1);
				}
			}
		}));
		riTimer.setCycleCount(Timeline.INDEFINITE);
		riTimer.play();
	}

	public static void main(String[] args) {
		Application.launch(args);
	}

}
