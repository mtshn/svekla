package ru.ac.phyche.gcms.svekla.javafxgui;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.concurrent.Worker;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
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
import ru.ac.phyche.gcms.ei2fp_java.Inference;

public class JavaFXGUI extends Application {
	public static final int MAX_ISOMERS_SURGE = 3000000;
	public static final int MAX_MW_SURGE = 155;
	public static final int MAX_ISOMERS_TEXT_AREA = 1000000; // Chars
	public static final int MAX_ISOMERS_SMART_SORTED = 10000; // Lines (isomers)

	private class Sortable implements Comparable<Sortable> {
		public String s;
		public float value;

		@Override
		public int compareTo(Sortable o) {
			return ((Float) value).compareTo(o.value);
		}

	}

	@Override
	public void start(Stage primaryStage) throws Exception {
		primaryStage.setTitle("SVEKLA: A GUI for GC-MS retention index and mass spectra prediction");
		FXMLLoader loader = new FXMLLoader();
		loader.setLocation(new File("./svekla.fxml").toURI().toURL());
		VBox vBox = loader.<VBox>load();
		Scene scene = new Scene(vBox, 1350, 900);
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
		TextField molecularFormulaInput = (TextField) vBox.lookup("#tf_molecular_formula_input");
		TextField tfSurgeParams = (TextField) vBox.lookup("#tf_surge_params");
		TextField tfWeightPolarRI = (TextField) vBox.lookup("#tf_weight_wax");
		TextField tfThresholdNonPolarRI = (TextField) vBox.lookup("#tf_non_polar_threshold");
		TextField tfThresholdPolarRI = (TextField) vBox.lookup("#tf_polar_threshold");

		TextArea targetMS = (TextArea) vBox.lookup("#textareaspectrum");
		TextArea textareaMSMS = (TextArea) vBox.lookup("#msmstextarea");
		TextArea textareaOutIsomers = (TextArea) vBox.lookup("#ta_output");

		CheckBox checkbox1 = (CheckBox) vBox.lookup("#checkboxsquareroot");
		CheckBox checkbox2 = (CheckBox) vBox.lookup("#checkboxallpeaks");
		CheckBox checkboxMSMS1 = (CheckBox) vBox.lookup("#checkboxmsmspos");
		CheckBox checkboxMSMS2 = (CheckBox) vBox.lookup("#checkboxmsmsneg");
		CheckBox checkboxMSMS3 = (CheckBox) vBox.lookup("#checkboxmsmsnewwin");
		CheckBox checkboxMSMS4 = (CheckBox) vBox.lookup("#checkboxmsmslargest");
		CheckBox cbNoTriple = (CheckBox) vBox.lookup("#cb_no_triple");
		CheckBox cbNoSmallRings = (CheckBox) vBox.lookup("#cb_no_small_rings");
		CheckBox cbNoX_X_X = (CheckBox) vBox.lookup("#cb_no_x_x_x");
		CheckBox cbCustomSurge = (CheckBox) vBox.lookup("#cb_custom_surge");

		WebView webView = (WebView) vBox.lookup("#webview");
		Canvas canvas1 = (Canvas) vBox.lookup("#canvas1");
		Canvas canvas2 = (Canvas) vBox.lookup("#canvas2");
		Button button = (Button) vBox.lookup("#button");
		Button button2 = (Button) vBox.lookup("#button2");
		Button buttonLoad = (Button) vBox.lookup("#loadbutton");
		Button buttonMSMS = (Button) vBox.lookup("#buttonmsms");
		Button buttonUseSpectrum = (Button) vBox.lookup("#usepastebutton");
		Button buttonCloseSpectra = (Button) vBox.lookup("#closespectra");
		Button buttonPubchemIsomers = (Button) vBox.lookup("#b_from_pubchem");
		Button buttonSurgeIsomers = (Button) vBox.lookup("#b_from_surge");
		Button buttonSmartFPSort = (Button) vBox.lookup("#b_sort_fp");

		Button buttonPredictRI = (Button) vBox.lookup("#b_predict_ri");
		Button buttonFilterRI = (Button) vBox.lookup("#b_filter_ri");
		Button buttonSortRI = (Button) vBox.lookup("#b_sort_ri");

		StringBuffer isomers = new StringBuffer();

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
		RIPrediction rp = new RIPrediction();

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

		buttonPubchemIsomers.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String url = molecularFormulaInput.getText().trim();
				url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/" + url;
				url = url + "/property/IsomericSMILES/TXT";
				String s = "";
				try {
					BufferedReader reader = new BufferedReader(new InputStreamReader(new URL(url).openStream()));
					String s1 = reader.readLine();
					while (s1 != null) {
						if (!s1.trim().equals("")) {
							if ((!s1.contains(".")) && (!s1.contains("[13C")) && (!s1.contains("[2H]"))) {
								s = s + s1 + "\n";
							}
						}
						s1 = reader.readLine();
					}
					reader.close();
				} catch (Exception e) {
					e.printStackTrace();
					Alert a = new Alert(AlertType.ERROR);
					String st = ExceptionUtils.getStackTrace(e);
					a.setContentText(st);
					a.show();
				}
				isomers.setLength(0);
				isomers.append(s);
				System.out.println("N isomers " + isomers.toString().split("\\n").length);
				System.out.println("Done");
				if (isomers.length() > MAX_ISOMERS_TEXT_AREA) {
					textareaOutIsomers.setText(isomers.substring(0, MAX_ISOMERS_TEXT_AREA)
							+ "\n....\n....too much isomers! Isomers are stored and can be sorted but cannot be shown");
				} else {
					textareaOutIsomers.setText(isomers.toString());
				}
			}
		});

		buttonSurgeIsomers.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String s = "";
				try {
					String molFormula = molecularFormulaInput.getText().trim();
					float mass = (float) MolecularFormulaManipulator.getMass(MolecularFormulaManipulator
							.getMolecularFormula(molFormula, DefaultChemObjectBuilder.getInstance()));
					if (mass < MAX_MW_SURGE) {
						boolean noTriple = cbNoTriple.selectedProperty().getValue();
						boolean noSmallRings = cbNoSmallRings.selectedProperty().getValue();
						boolean noX_X_X = cbNoX_X_X.selectedProperty().getValue();
						boolean customSurge = cbCustomSurge.selectedProperty().getValue();
						String customSurgeLine = tfSurgeParams.getText();

						int nIsomers = SurgeIsomers.surgeCount(molFormula, noTriple, noSmallRings, noX_X_X, customSurge,
								customSurgeLine);
						System.out.println(nIsomers);
						if (nIsomers < MAX_ISOMERS_SURGE) {
							String[] isomers = SurgeIsomers.surgeGenerate(molFormula, noTriple, noSmallRings, noX_X_X,
									customSurge, customSurgeLine);
							s = SurgeIsomers.concatNewLine(isomers);
						} else {
							Alert a = new Alert(AlertType.ERROR);
							String st = "Surge generated " + nIsomers
									+ " for this molecular formula! Too many isomers!";
							a.setContentText(st);
							a.setHeaderText("Too many isomers!");
							a.show();
						}
					} else {
						Alert a = new Alert(AlertType.ERROR);
						String st = "Surge can be used only with molecular weight < " + MAX_MW_SURGE;
						a.setContentText(st);
						a.setHeaderText("Too many isomers!");
						a.show();
					}
				} catch (Exception e) {
					e.printStackTrace();
					Alert a = new Alert(AlertType.ERROR);
					String st = ExceptionUtils.getStackTrace(e);
					a.setContentText(st);
					a.show();
				}
				isomers.setLength(0);
				isomers.append(s);
				System.out.println("N isomers " + isomers.toString().split("\\n").length);
				System.out.println("Done");
				if (isomers.length() > MAX_ISOMERS_TEXT_AREA) {
					textareaOutIsomers.setText(isomers.substring(0, MAX_ISOMERS_TEXT_AREA)
							+ "\n....\n....too much isomers! Isomers are stored and can be sorted but cannot be shown");
				} else {
					textareaOutIsomers.setText(isomers.toString());
				}
			}
		});

		buttonSmartFPSort.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String[] splt = isomers.toString().split("\\n");
				String[] smiles = new String[splt.length];
				for (int i = 0; i < splt.length; i++) {
					String t = splt[i].trim();
					if (!t.equals("")) {
						smiles[i] = t.split("\\s+")[0];
					}
				}
				float[] fp0 = MSpectra.fingerprintsFromSpectrum(target[0]);
				float[] crossEntropies = Inference.crossentropy(fp0, smiles);
				Sortable[] s3 = new Sortable[splt.length];
				for (int i = 0; i < splt.length; i++) {
					String t = splt[i].trim();
					if (!t.equals("")) {
						splt[i] = splt[i].trim() + " S: " + crossEntropies[i];
						s3[i] = new Sortable();
						s3[i].s = splt[i];
						s3[i].value = crossEntropies[i];
					} else {
						s3[i] = new Sortable();
						s3[i].s = "";
						s3[i].value = Float.POSITIVE_INFINITY;
					}

				}
				Arrays.sort(s3);
				String result = "";
				for (int i = 0; (i < splt.length) && (i < MAX_ISOMERS_SMART_SORTED); i++) {
					result = result + s3[i].s + "\n";
				}
				textareaOutIsomers.setText(result);
				isomers.setLength(0);
				isomers.append(result);
			}

		});

		buttonPredictRI.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String[] splt = isomers.toString().split("\\n");
				String result = "";
				if (splt.length < MAX_ISOMERS_SMART_SORTED + 10) {
					for (int i = 0; i < splt.length; i++) {
						String t = splt[i].trim();
						if (!t.equals("")) {
							t = t.split("\\s+")[0];
							float[] ri = rp.predictNonpolarPolar(t);
							result = result + splt[i] + " NP: " + ri[0] + " WAX: " + ri[1] + "\n";
						} else {
							result = result + "\n";
						}

					}
				} else {
					Alert a = new Alert(AlertType.ERROR);
					String st = "The maximum number of structures for which retention indices can be predicted is "
							+ MAX_ISOMERS_SMART_SORTED;
					a.setContentText(st);
					a.setHeaderText("Too many structures!");
					a.show();
				}
				textareaOutIsomers.setText(result);
				isomers.setLength(0);
				isomers.append(result);
			}

		});

		buttonSortRI.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String[] splt = isomers.toString().split("\\n");
				float targetPolar = Float.parseFloat(polarTarget.getText().trim());
				float targetNonPolar = Float.parseFloat(nonpolarTarget.getText().trim());
				float weightPolar = Float.parseFloat(tfWeightPolarRI.getText().trim());
				Sortable[] sort = new Sortable[splt.length];
				if (splt.length < MAX_ISOMERS_SMART_SORTED + 10) {
					for (int i = 0; i < splt.length; i++) {
						String t = splt[i].trim();
						float value = Float.POSITIVE_INFINITY;
						if (!t.equals("")) {
							String[] x = t.split("\\s+");
							float nonpolar = -1000;
							float polar = -1000;
							for (int j = 0; j < x.length; j++) {
								if (x[j].equals("NP:")) {
									nonpolar = Float.parseFloat(x[j + 1]);
								}
								if (x[j].equals("WAX:")) {
									polar = Float.parseFloat(x[j + 1]);
								}
							}
							value = Math.abs(nonpolar - targetNonPolar) + weightPolar * Math.abs(polar - targetPolar);
							if ((nonpolar < -900) || (polar < -900)) {
								Alert a = new Alert(AlertType.ERROR);
								String st = "Predict retention indices before sorting ";
								a.setContentText(st);
								a.show();
								throw new RuntimeException("Predict retention indices before sorting");
							}
						}
						sort[i] = new Sortable();
						sort[i].value = value;
						sort[i].s = t;
					}
				} else {
					Alert a = new Alert(AlertType.ERROR);
					String st = "The maximum number of structures for which retention indices can be predicted is "
							+ MAX_ISOMERS_SMART_SORTED;
					a.setContentText(st);
					a.setHeaderText("Too many structures!");
					a.show();
				}
				Arrays.sort(sort);
				String result = "";
				for (int i = 0; (i < splt.length) && (i < MAX_ISOMERS_SMART_SORTED); i++) {
					result = result + sort[i].s + " D_RI " + sort[i].value + "\n";
				}
				textareaOutIsomers.setText(result);
				isomers.setLength(0);
				isomers.append(result);
			}
		});

		buttonFilterRI.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent actionEvent) {
				String[] splt = isomers.toString().split("\\n");
				float targetPolar = Float.parseFloat(polarTarget.getText().trim());
				float targetNonPolar = Float.parseFloat(nonpolarTarget.getText().trim());
				float thresholdNonPolar = Float.parseFloat(tfThresholdNonPolarRI.getText().trim());
				float thresholdPolar = Float.parseFloat(tfThresholdPolarRI.getText().trim());

				String result = "";
				if (splt.length < MAX_ISOMERS_SMART_SORTED + 10) {
					for (int i = 0; i < splt.length; i++) {
						String t = splt[i].trim();
						if (!t.equals("")) {
							String[] x = t.split("\\s+");
							float nonpolar = -1000;
							float polar = -1000;
							for (int j = 0; j < x.length; j++) {
								if (x[j].equals("NP:")) {
									nonpolar = Float.parseFloat(x[j + 1]);
								}
								if (x[j].equals("WAX:")) {
									polar = Float.parseFloat(x[j + 1]);
								}
							}
							if ((nonpolar < -900) || (polar < -900)) {
								Alert a = new Alert(AlertType.ERROR);
								String st = "Predict retention indices before sorting ";
								a.setContentText(st);
								a.show();
								throw new RuntimeException("Predict retention indices before sorting");
							}
							boolean o = true;
							if (Math.abs(nonpolar - targetNonPolar) > thresholdNonPolar) {
								o = false;
							}
							if (Math.abs(polar - targetPolar) > thresholdPolar) {
								o = false;
							}
							if (o) {
								result = result + t + "\n";
							}
						}
					}
				} else {
					Alert a = new Alert(AlertType.ERROR);
					String st = "The maximum number of structures for which retention indices can be predicted is "
							+ MAX_ISOMERS_SMART_SORTED;
					a.setContentText(st);
					a.setHeaderText("Too many structures!");
					a.show();
				}
				textareaOutIsomers.setText(result);
				isomers.setLength(0);
				isomers.append(result);
			}
		});

		primaryStage.show();

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
